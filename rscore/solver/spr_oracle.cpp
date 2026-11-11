// Restricted SPR clone oracle for the reliable-quartet solver.
//
// This intentionally does not change the solver search path.  It applies one
// non-root-child SPR, rebuilds the clone bitsets/sizes, and verifies that
// recomputing only the theoretically affected cache rows gives exactly the
// same delta as a full clone refresh.
#define SOLVER_NO_MAIN
#include "solver.cpp"

#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>

namespace {

struct SPRMove {
    int cutParent = -1; // p
    int cutChild = -1;  // r, the subtree being moved
    int targetParent = -1; // q
    int targetChild = -1;  // t, old child of q
};

vector<GeneTree> readGenes(string const& path) {
    std::ifstream fin(path);
    if (!fin) { std::cerr << "Error: cannot open gene file " << path << '\n'; std::exit(2); }
    vector<GeneTree> genes;
    string line;
    while (std::getline(fin, line)) {
        if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
        newick::Newick nw(line);
        genes.push_back(buildGene(nw.root()));
    }
    if (genes.empty()) { std::cerr << "Error: no gene trees\n"; std::exit(2); }
    return genes;
}

std::unordered_set<string> readTaxonNames(string const& path) {
    std::ifstream fin(path);
    if (!fin) { std::cerr << "Error: cannot open gene file " << path << '\n'; std::exit(2); }
    std::unordered_set<string> names;
    string line;
    while (std::getline(fin, line)) {
        if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
        newick::Newick nw(line);
        vector<newick::Newick::Node const*> st{nw.root()};
        while (!st.empty()) {
            auto const* u = st.back(); st.pop_back();
            if (u->children.empty()) names.insert(u->name);
            else for (auto const* c : u->children) st.push_back(c);
        }
    }
    return names;
}

SpeciesTree readSpecies(string const& path, size_t nTaxa,
                        std::unordered_set<string> const& expectedNames) {
    std::ifstream fin(path);
    if (!fin) { std::cerr << "Error: cannot open species file " << path << '\n'; std::exit(2); }
    string line;
    while (std::getline(fin, line)) {
        if (line.find_first_not_of(" \t\r\n") != string::npos) {
            newick::Newick nw(line);
            std::unordered_set<string> speciesNames;
            vector<newick::Newick::Node const*> st{nw.root()};
            while (!st.empty()) {
                auto const* u = st.back(); st.pop_back();
                if (u->children.empty()) speciesNames.insert(u->name);
                else for (auto const* c : u->children) st.push_back(c);
            }
            if (speciesNames != expectedNames) {
                std::cerr << "Error: species/gene taxon sets differ\n";
                std::exit(2);
            }
            SpeciesTree tree = buildSpecies(nw.root(), nTaxa);
            return tree;
        }
    }
    std::cerr << "Error: no species tree\n";
    std::exit(2);
}

void rebuildBits(SpeciesTree& t) {
    std::fill(t.bits.begin(), t.bits.end(), 0);
    std::fill(t.sz.begin(), t.sz.end(), 0);
    std::function<void(int)> visit = [&](int u) {
        if (t.taxon[u] >= 0) {
            t.bits[(size_t)u * t.W + ((size_t)t.taxon[u] >> 6)] |= 1ULL << (t.taxon[u] & 63);
            t.sz[u] = 1;
            return;
        }
        if (t.c0[u] < 0 || t.c1[u] < 0) {
            std::cerr << "internal SPR error: disconnected node " << u << '\n';
            std::exit(1);
        }
        visit(t.c0[u]);
        visit(t.c1[u]);
        u64* dst = t.bits.data() + (size_t)u * t.W;
        u64 const* a = t.bits.data() + (size_t)t.c0[u] * t.W;
        u64 const* b = t.bits.data() + (size_t)t.c1[u] * t.W;
        for (size_t w = 0; w < t.W; ++w) dst[w] = a[w] | b[w];
        t.sz[u] = t.sz[t.c0[u]] + t.sz[t.c1[u]];
    };
    visit(t.root);
}

vector<i64> refreshCache(SpeciesTree const& t, vector<GeneTree> const& genes,
                         int T, int threads) {
    int maxN = 1;
    for (auto const& g : genes) maxN = std::max(maxN, g.n);
    vector<Scratch> scratch(std::max(1, threads));
    for (auto& sc : scratch) { sc.resize(maxN); sc.t1.assign(t.W, 0); sc.t2.assign(t.W, 0); }
    vector<i64> cache((size_t)t.internal.size() * T, 0);
    #pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int k = 0; k < (int)t.internal.size(); ++k) {
        Scratch& sc = scratch[omp_get_thread_num()];
        int const v = t.internal[(size_t)k];
        int const a = t.c0[v], b = t.c1[v];
        u64 const* ba = t.bits.data() + (size_t)a * t.W;
        u64 const* bb = t.bits.data() + (size_t)b * t.W;
        for (size_t w = 0; w < t.W; ++w) sc.t1[w] = ~(ba[w] | bb[w]);
        for (int g = 0; g < T; ++g)
            cache[(size_t)k * T + g] = Wnode(genes[(size_t)g], ba, bb, sc.t1.data(),
                t.sz[a], t.sz[b], (int)t.sz[t.root] - t.sz[a] - t.sz[b], t.W, sc);
    }
    return cache;
}

bool isDescendant(SpeciesTree const& t, int node, int ancestor) {
    for (int u = node; u >= 0; u = t.par[u]) if (u == ancestor) return true;
    return false;
}

bool validMove(SpeciesTree const& t, SPRMove const& m) {
    if (m.cutParent < 0 || m.cutChild < 0 || m.targetParent < 0 || m.targetChild < 0) return false;
    if (m.cutParent >= t.n || m.cutChild >= t.n || m.targetParent >= t.n || m.targetChild >= t.n) return false;
    if (t.par[m.cutChild] != m.cutParent || t.par[m.targetChild] != m.targetParent) return false;
    if (t.c0[m.cutParent] < 0 || t.c0[m.targetParent] < 0) return false;
    int const gp = t.par[m.cutParent];
    if (gp < 0 || (t.c0[gp] != m.cutParent && t.c1[gp] != m.cutParent)) return false;
    // The target edge must lie outside R.  Check both endpoints as subsets of R.
    if (isDescendant(t, m.targetParent, m.cutChild)
        || isDescendant(t, m.targetChild, m.cutChild)) return false;
    if (m.targetChild == m.cutParent) return false; // edge removed with p

    int const sibling = (t.c0[m.cutParent] == m.cutChild) ? t.c1[m.cutParent] : t.c0[m.cutParent];
    if (m.targetParent == gp && m.targetChild == sibling) return false; // no-op
    if (m.targetParent == m.cutParent || m.targetChild == m.cutParent) return false;
    return true;
}

vector<char> affectedMask(SpeciesTree const& t, SPRMove const& m) {
    vector<char> affected((size_t)t.n, 0);
    affected[(size_t)m.cutParent] = 1; // p is repurposed connector
    for (int u = t.par[m.cutParent]; u >= 0; u = t.par[u]) affected[(size_t)u] = 1;
    for (int u = m.targetParent; u >= 0; u = t.par[u]) affected[(size_t)u] = 1;
    return affected;
}

void applySPR(SpeciesTree& t, SPRMove const& m) {
    if (!validMove(t, m)) { std::cerr << "invalid SPR move\n"; std::exit(1); }
    int const p = m.cutParent, r = m.cutChild;
    int const gp = t.par[p];
    int const s = (t.c0[p] == r) ? t.c1[p] : t.c0[p];
    // Suppress p in the old location.
    if (t.c0[gp] == p) t.c0[gp] = s;
    else if (t.c1[gp] == p) t.c1[gp] = s;
    else { std::cerr << "SPR parent link mismatch\n"; std::exit(1); }
    t.par[s] = gp;
    // Reuse p as the connector on q -> t.
    int const q = m.targetParent, target = m.targetChild;
    if (t.c0[q] == target) t.c0[q] = p;
    else if (t.c1[q] == target) t.c1[q] = p;
    else { std::cerr << "SPR target link mismatch\n"; std::exit(1); }
    t.par[p] = q;
    t.c0[p] = target; t.c1[p] = r;
    t.par[target] = p; t.par[r] = p;
    rebuildBits(t);
}

vector<i64> refreshCacheRows(SpeciesTree const& t, vector<GeneTree> const& genes,
                             int T, int threads, vector<char> const& affected,
                             vector<i64> const& oldCache) {
    int maxN = 1;
    for (auto const& g : genes) maxN = std::max(maxN, g.n);
    vector<Scratch> scratch(std::max(1, threads));
    for (auto& sc : scratch) { sc.resize(maxN); sc.t1.assign(t.W, 0); sc.t2.assign(t.W, 0); }
    vector<i64> cache = oldCache;
    #pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int k = 0; k < (int)t.internal.size(); ++k) {
        int const v = t.internal[(size_t)k];
        if (!affected[(size_t)v]) continue;
        Scratch& sc = scratch[omp_get_thread_num()];
        int const a = t.c0[v], b = t.c1[v];
        u64 const* ba = t.bits.data() + (size_t)a * t.W;
        u64 const* bb = t.bits.data() + (size_t)b * t.W;
        for (size_t w = 0; w < t.W; ++w) sc.t1[w] = ~(ba[w] | bb[w]);
        for (int g = 0; g < T; ++g)
            cache[(size_t)k * T + g] = Wnode(genes[(size_t)g], ba, bb, sc.t1.data(),
                t.sz[a], t.sz[b], (int)t.sz[t.root] - t.sz[a] - t.sz[b], t.W, sc);
    }
    return cache;
}

bool checkInvariants(SpeciesTree const& t, int nTaxa) {
    if (t.root < 0 || t.root >= t.n || t.par[t.root] != -1) return false;
    vector<char> seen((size_t)t.n, 0), taxaSeen((size_t)nTaxa, 0);
    vector<int> leaves;
    std::function<void(int,int)> visit = [&](int u, int parent) {
        if (u < 0 || u >= t.n || seen[(size_t)u]) throw 1;
        seen[(size_t)u] = 1;
        if (t.par[u] != parent) throw 1;
        if (t.taxon[u] >= 0) {
            if (t.taxon[u] >= nTaxa || taxaSeen[(size_t)t.taxon[u]]) throw 1;
            if (t.c0[u] >= 0 || t.c1[u] >= 0 || t.sz[u] != 1) throw 1;
            for (size_t w = 0; w < t.W; ++w) {
                u64 const want = (w == ((size_t)t.taxon[u] >> 6)) ? (1ULL << (t.taxon[u] & 63)) : 0;
                if (t.bits[(size_t)u * t.W + w] != want) throw 1;
            }
            taxaSeen[(size_t)t.taxon[u]] = 1;
            leaves.push_back(t.taxon[u]);
            return;
        }
        if (t.c0[u] < 0 || t.c1[u] < 0 || t.c0[u] == t.c1[u]) throw 1;
        visit(t.c0[u], u); visit(t.c1[u], u);
        if (t.sz[u] != t.sz[t.c0[u]] + t.sz[t.c1[u]]) throw 1;
        for (size_t w = 0; w < t.W; ++w) {
            u64 const want = t.bits[(size_t)t.c0[u] * t.W + w] | t.bits[(size_t)t.c1[u] * t.W + w];
            if (t.bits[(size_t)u * t.W + w] != want) throw 1;
        }
    };
    try { visit(t.root, -1); } catch (...) { return false; }
    for (char x : seen) if (!x) return false;
    for (char x : taxaSeen) if (!x) return false;
    return (int)leaves.size() == nTaxa && t.sz[t.root] == nTaxa
        && vector<int>(t.internal) == [&] { vector<int> v; for (int u = 0; u < t.n; ++u) if (t.c0[u] >= 0) v.push_back(u); return v; }();
}

__int128 sumCache(vector<i64> const& c) {
    __int128 total = 0;
    for (i64 x : c) total += x;
    return total;
}

string printInt128(__int128 x) {
    if (x == 0) return "0";
    bool neg = x < 0;
    if (neg) x = -x;
    string s;
    while (x > 0) { s.push_back(char('0' + x % 10)); x /= 10; }
    if (neg) s.push_back('-');
    std::reverse(s.begin(), s.end());
    return s;
}

struct BatchResult {
    __int128 detachDelta = 0;
    i64 base = 0;
    vector<i64> dA, dB, dC, prefix, newNode, detachedW;
    vector<char> main;
};

void markMain(SpeciesTree const& t, int u, vector<char>& main) {
    if (u < 0 || main[(size_t)u]) return;
    main[(size_t)u] = 1;
    if (t.c0[u] >= 0) {
        markMain(t, t.c0[u], main);
        markMain(t, t.c1[u], main);
    }
}

SpeciesTree detachedTree(SpeciesTree const& base, SPRMove const& m) {
    SpeciesTree d = base;
    int const p = m.cutParent, r = m.cutChild, gp = d.par[p];
    int const s = (d.c0[p] == r) ? d.c1[p] : d.c0[p];
    if (d.c0[gp] == p) d.c0[gp] = s;
    else d.c1[gp] = s;
    d.par[s] = gp;
    d.par[p] = -1;
    d.c0[p] = d.c1[p] = -1;
    d.sz[p] = 0;
    std::fill(d.bits.begin() + (size_t)p * d.W,
              d.bits.begin() + (size_t)(p + 1) * d.W, 0);
    return d;
}

void rebuildMainBits(SpeciesTree& t, vector<char> const& main) {
    std::function<void(int)> visit = [&](int u) {
        if (!main[(size_t)u]) return;
        if (t.taxon[u] >= 0) {
            std::fill(t.bits.begin() + (size_t)u * t.W,
                      t.bits.begin() + (size_t)(u + 1) * t.W, 0);
            t.bits[(size_t)u * t.W + ((size_t)t.taxon[u] >> 6)] |= 1ULL << (t.taxon[u] & 63);
            t.sz[u] = 1;
            return;
        }
        visit(t.c0[u]); visit(t.c1[u]);
        u64* dst = t.bits.data() + (size_t)u * t.W;
        u64 const* a = t.bits.data() + (size_t)t.c0[u] * t.W;
        u64 const* b = t.bits.data() + (size_t)t.c1[u] * t.W;
        for (size_t w = 0; w < t.W; ++w) dst[w] = a[w] | b[w];
        t.sz[u] = t.sz[t.c0[u]] + t.sz[t.c1[u]];
    };
    visit(t.root);
}

BatchResult generalizedBatch(SpeciesTree const& base, SPRMove const& m,
                             vector<GeneTree> const& genes, int T, int threads,
                             vector<i64> const& oldCache, int nTaxa) {
    SpeciesTree d = detachedTree(base, m);
    BatchResult out;
    out.main.assign((size_t)d.n, 0);
    markMain(d, d.root, out.main);
    rebuildMainBits(d, out.main);
    size_t const W = d.W;
    int const r = m.cutChild;
    int const rsize = base.sz[r];
    vector<u64> block(base.bits.begin() + (size_t)r * W,
                      base.bits.begin() + (size_t)(r + 1) * W);
    vector<u64> full(W, ~0ULL);
    if (nTaxa % 64) full.back() = (1ULL << (nTaxa % 64)) - 1;
    out.dA.assign((size_t)d.n, 0); out.dB.assign((size_t)d.n, 0);
    out.dC.assign((size_t)d.n, 0); out.prefix.assign((size_t)d.n, 0);
    out.newNode.assign((size_t)d.n, 0); out.detachedW.assign((size_t)d.n, 0);
    int maxN = 1;
    for (auto const& g : genes) maxN = std::max(maxN, g.n);
    vector<Scratch> scratch(std::max(1, threads));
    for (auto& sc : scratch) { sc.resize(maxN); sc.t1.assign(W, 0); sc.t2.assign(W, 0); }
    #pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int v = 0; v < d.n; ++v) {
        if (!out.main[(size_t)v] || d.c0[v] < 0) continue;
        Scratch& sc = scratch[omp_get_thread_num()];
        int const a = d.c0[v], b = d.c1[v];
        u64 const* ba = d.bits.data() + (size_t)a * W;
        u64 const* bb = d.bits.data() + (size_t)b * W;
        vector<u64> qa(W), qb(W), qc(W);
        for (size_t w = 0; w < W; ++w) { qa[w] = ba[w]; qb[w] = bb[w]; qc[w] = full[w] & ~(qa[w] | qb[w]); }
        i64 oldW = 0;
        for (int g = 0; g < T; ++g)
            oldW += Wnode(genes[(size_t)g], qa.data(), qb.data(), qc.data(),
                          d.sz[a], d.sz[b], nTaxa - d.sz[a] - d.sz[b], W, sc);
        out.detachedW[(size_t)v] = oldW;
        for (size_t w = 0; w < W; ++w) { qa[w] = ba[w] | block[w]; qb[w] = bb[w]; qc[w] = full[w] & ~(qa[w] | qb[w]); }
        i64 wA = 0;
        for (int g = 0; g < T; ++g)
            wA += Wnode(genes[(size_t)g], qa.data(), qb.data(), qc.data(),
                        d.sz[a] + rsize, d.sz[b], nTaxa - d.sz[a] - rsize - d.sz[b], W, sc);
        out.dA[(size_t)v] = wA - oldW;
        for (size_t w = 0; w < W; ++w) { qa[w] = ba[w]; qb[w] = bb[w] | block[w]; qc[w] = full[w] & ~(qa[w] | qb[w]); }
        i64 wB = 0;
        for (int g = 0; g < T; ++g)
            wB += Wnode(genes[(size_t)g], qa.data(), qb.data(), qc.data(),
                        d.sz[a], d.sz[b] + rsize, nTaxa - d.sz[a] - d.sz[b] - rsize, W, sc);
        out.dB[(size_t)v] = wB - oldW;
        for (size_t w = 0; w < W; ++w) { qa[w] = ba[w]; qb[w] = bb[w]; qc[w] = full[w] & ~(qa[w] | qb[w]); }
        i64 wC = 0;
        for (int g = 0; g < T; ++g)
            wC += Wnode(genes[(size_t)g], qa.data(), qb.data(), qc.data(),
                        d.sz[a], d.sz[b], nTaxa - d.sz[a] - d.sz[b], W, sc);
        out.dC[(size_t)v] = wC - oldW;
    }
    vector<int> idx(base.n, -1);
    for (int k = 0; k < (int)base.internal.size(); ++k) idx[(size_t)base.internal[(size_t)k]] = k;
    for (int v = 0; v < d.n; ++v) if (out.main[(size_t)v] && d.c0[v] >= 0) {
        out.base += out.dC[(size_t)v];
        int const k = idx[(size_t)v];
        if (k >= 0) {
            __int128 oldAgg = 0;
            for (int g = 0; g < T; ++g) oldAgg += oldCache[(size_t)k * T + g];
            out.detachDelta += (__int128)out.detachedW[(size_t)v] - oldAgg;
        }
    }
    // Remove the old p node, whose cache row is not part of the detached main tree.
    int const pk = idx[(size_t)m.cutParent];
    for (int g = 0; g < T; ++g) out.detachDelta -= oldCache[(size_t)pk * T + g];
    vector<int> st{d.root};
    while (!st.empty()) {
        int const v = st.back(); st.pop_back();
        if (!out.main[(size_t)v] || d.c0[v] < 0) continue;
        int const a = d.c0[v], b = d.c1[v];
        out.prefix[(size_t)a] = out.prefix[(size_t)v] + out.dA[(size_t)v] - out.dC[(size_t)v];
        out.prefix[(size_t)b] = out.prefix[(size_t)v] + out.dB[(size_t)v] - out.dC[(size_t)v];
        st.push_back(a); st.push_back(b);
    }
    #pragma omp parallel for schedule(dynamic) num_threads(threads)
    for (int t = 0; t < d.n; ++t) {
        if (!out.main[(size_t)t] || d.par[t] < 0) continue;
        Scratch& sc = scratch[omp_get_thread_num()];
        u64 const* bt = d.bits.data() + (size_t)t * W;
        vector<u64> qc(W);
        for (size_t w = 0; w < W; ++w) qc[w] = full[w] & ~(bt[w] | block[w]);
        i64 nw = 0;
        for (int g = 0; g < T; ++g)
            nw += Wnode(genes[(size_t)g], bt, block.data(), qc.data(),
                        d.sz[t], rsize, nTaxa - d.sz[t] - rsize, W, sc);
        out.newNode[(size_t)t] = nw;
    }
    return out;
}

__int128 batchDelta(BatchResult const& b, int target) {
    return b.detachDelta + b.base + b.prefix[(size_t)target] + b.newNode[(size_t)target];
}

vector<SPRMove> enumerateMoves(SpeciesTree const& t, int maxCuts, int targetsPerCut) {
    vector<SPRMove> out;
    int cuts = 0;
    for (int p : t.internal) {
        if (t.par[p] < 0) continue;
        int const children[2] = {t.c0[p], t.c1[p]};
        for (int r : children) {
            if (cuts++ >= maxCuts && maxCuts > 0) return out;
            int added = 0;
            for (int q : t.internal) {
                int const ts[2] = {t.c0[q], t.c1[q]};
                for (int target : ts) {
                    SPRMove m{p, r, q, target};
                    if (!validMove(t, m)) continue;
                    out.push_back(m);
                    if (targetsPerCut > 0 && ++added >= targetsPerCut) break;
                }
                if (targetsPerCut > 0 && added >= targetsPerCut) break;
            }
        }
    }
    return out;
}

} // namespace

#ifndef SPR_ORACLE_NO_MAIN
int main(int argc, char** argv) {
    string speciesFile, genesFile;
    int threads = 4, T = -1, maxCuts = 12, targetsPerCut = 32;
    for (int i = 1; i < argc; ++i) {
        string s = argv[i];
        auto next = [&]() -> string { if (++i >= argc) std::exit(2); return argv[i]; };
        if (s == "-s") speciesFile = next();
        else if (s == "-i") genesFile = next();
        else if (s == "-t") threads = std::stoi(next());
        else if (s == "--genes") T = std::stoi(next());
        else if (s == "--max-cuts") maxCuts = std::stoi(next());
        else if (s == "--targets-per-cut") targetsPerCut = std::stoi(next());
        else { std::cerr << "unknown arg " << s << '\n'; return 2; }
    }
    if (speciesFile.empty() || genesFile.empty()) {
        std::cerr << "usage: spr_oracle -s species.tree -i genes.trees [-t N] [--genes K]"
                     " [--max-cuts K] [--targets-per-cut K]\n";
        return 2;
    }
    if (threads < 1) {
        std::cerr << "Error: -t must be positive\n";
        return 2;
    }
    std::unordered_set<string> expectedNames = readTaxonNames(genesFile);
    vector<GeneTree> genes = readGenes(genesFile);
    size_t const nTaxa = common::taxonName2ID.nTaxa();
    SpeciesTree base = readSpecies(speciesFile, nTaxa, expectedNames);
    T = T < 0 ? (int)genes.size() : std::min<int>(T, genes.size());
    if (T < 1 || !checkInvariants(base, (int)nTaxa)) {
        std::cerr << "FAIL base species-tree invariants\n";
        return 1;
    }
    vector<i64> oldCache = refreshCache(base, genes, T, threads);
    auto moves = enumerateMoves(base, maxCuts, targetsPerCut);
    int passed = 0, failed = 0;
    __int128 minDelta = 0, maxDelta = 0;
    bool haveDelta = false;
    for (SPRMove const& m : moves) {
        vector<char> affected = affectedMask(base, m);
        SpeciesTree clone = base;
        applySPR(clone, m);
        bool moveOk = checkInvariants(clone, (int)nTaxa) && clone.internal == base.internal;
        if (!moveOk) {
            std::cerr << "FAIL clone species-tree invariants or internal IDs\n";
            ++failed;
            continue;
        }
        vector<i64> newCache = refreshCache(clone, genes, T, threads);
        vector<i64> localCache = refreshCacheRows(clone, genes, T, threads, affected, oldCache);
        for (size_t ci = 0; ci < newCache.size() && moveOk; ++ci) {
            int const k = (int)(ci / (size_t)T);
            int const v = base.internal[(size_t)k];
            bool const changed = affected[(size_t)v] != 0;
            if (localCache[ci] != newCache[ci] && changed) {
                std::cerr << "FAIL affected local row node=" << v << " gene=" << (ci % T) << '\n';
                moveOk = false;
            }
            if (!changed && localCache[ci] != oldCache[ci]) {
                std::cerr << "FAIL unaffected local row node=" << v << " gene=" << (ci % T) << '\n';
                moveOk = false;
            }
        }
        __int128 localDelta128 = 0;
        if (moveOk) {
            for (int k = 0; k < (int)base.internal.size(); ++k) {
                int const v = base.internal[(size_t)k];
                bool const changed = affected[(size_t)v] != 0;
                for (int g = 0; g < T; ++g) {
                    i64 const diff = newCache[(size_t)k * T + g] - oldCache[(size_t)k * T + g];
                    if (!changed && diff != 0) {
                        std::cerr << "FAIL unaffected cache node=" << v << " gene=" << g << '\n';
                        moveOk = false;
                        break;
                    }
                    if (changed) localDelta128 += diff;
                }
                if (!moveOk) break;
            }
        }
        if (!moveOk) {
            ++failed;
            continue;
        }
        __int128 const fullDelta = sumCache(newCache) - sumCache(oldCache);
        BatchResult const batch = generalizedBatch(base, m, genes, T, threads, oldCache, (int)nTaxa);
        __int128 const batchDeltaValue = batchDelta(batch, m.targetChild);
        if (batchDeltaValue != fullDelta || (batchDeltaValue & 1) != 0) {
            std::cerr << "FAIL batch delta=" << printInt128(batchDeltaValue)
                      << " full=" << printInt128(fullDelta) << '\n';
            ++failed;
            continue;
        }
        if (localDelta128 != fullDelta || (localDelta128 & 1) != 0) {
            std::cerr << "FAIL delta local=" << printInt128(localDelta128)
                      << " full=" << printInt128(fullDelta) << '\n';
            ++failed;
            continue;
        }
        ++passed;
        if (!haveDelta) { minDelta = maxDelta = localDelta128; haveDelta = true; }
        else { minDelta = std::min(minDelta, localDelta128); maxDelta = std::max(maxDelta, localDelta128); }
    }
    std::cout << "SPR oracle: moves=" << moves.size() << " passed=" << passed
              << " failed=" << failed;
    if (haveDelta) std::cout << " delta_range=[" << printInt128(minDelta) << ',' << printInt128(maxDelta) << ']';
    std::cout << " genes=" << T << "\n";
    return failed ? 1 : 0;
}
#endif // SPR_ORACLE_NO_MAIN
