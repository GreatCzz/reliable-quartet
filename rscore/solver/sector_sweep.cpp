// sector_sweep.cpp — phase 1 of the exact sector ("leaf placement surface") sweep,
// plus a parallel-scaling benchmark.
//
// For every internal node m of the species tree S and every leaf x this computes
//
//     T_x(m,d) = # reliable quartets  x a | b c  whose partner a lies in direction d of m
//
// via Theorem 1 of rscore/SECTOR_MATH.md:
//
//     T_g(m,d)   = SUM_q  P_d(q,i_q(x)) * CPS_d(q,i_q(x))
//     P_d(q,i)   = |D_d ∩ C_i(q) ∩ A_x|
//     CPS_d(q,i) = SUM_{j<l, j,l≠i} ( Q_d'(q,j)Q_d''(q,l) + Q_d'(q,l)Q_d''(q,j) )
//     Q_d(q,j)   = |D_d ∩ C_j(q) ∩ A_q|
//
// Theorem 3 makes the whole x-array cost one pass per (m, gene): P factorises as
// |D_d ∩ C_i| - |gamma^c ∩ D_d ∩ C_i|, so it depends on x only through (C_i, gamma),
// and the per-x accumulation is the edge-difference path sum
//
//     T_x(m,d) = TOTALP_d + SUM over path edges (p->u) of [ EDG_d(p->u) - GP_d(p) ]
//
// Correction required because the tripartition of S-m still contains x: the true one,
// for (S-x)-m, is D_c \ {x} in the direction c containing x, hence
//
//     T_x(m,c) -= S_x(m,c),  S_x(m,d) = SUM_q CPS_d(q,i_q(x))   (same sweep, P := 1)
//
// Strategies (benchmark):
//   m       task = one m,                 parallel over m
//   block   task = B consecutive m's,     gene-major inner loop, parallel over blocks
//   nested  all threads share one m,      parallel over genes, per-thread accumulators
//
// Build: g++ -std=c++20 -O3 -march=native -fopenmp sector_sweep.cpp -o sector_sweep

#include <algorithm>
#include <chrono>
#include <functional>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include <omp.h>

// ---------------------------------------------------------------- newick ----
struct Gene {
    int n = 0;
    std::vector<int> head, ch, par, tax, post, depth, clId, slot;
    std::vector<int> clHead, clLeaf;          // CSR: leaves of each solid cluster
    std::vector<char> solid, small, inGamma;
    int nleaf = 0, maxdeg = 0, ncl = 0, bestCl = 0, nslot = 0;
};

static void parse_newick(const std::string &s, Gene &g)
{
    std::vector<std::vector<int>> kids;
    std::vector<int> taxv, solidv, stk;
    int root = -1;
    size_t i = 0, n = s.size();
    while (i < n) {
        char c = s[i];
        if (c == '(') {
            kids.emplace_back(); taxv.push_back(-1); solidv.push_back(0);
            stk.push_back((int)kids.size() - 1); ++i;
        } else if (c == ',' || c == ';' || c == ' ' || c == '\t' || c == '\r' || c == '\n') {
            ++i;
        } else if (c == ')') {
            int u = stk.back(); stk.pop_back(); ++i;
            if (i < n && (s[i] == 'S' || s[i] == 'D')) { solidv[u] = (s[i] == 'S'); ++i; }
            if (stk.empty()) root = u; else kids[stk.back()].push_back(u);
        } else if (c >= '0' && c <= '9') {
            size_t j = i; int t = 0;
            while (j < n && s[j] >= '0' && s[j] <= '9') { t = t * 10 + (s[j] - '0'); ++j; }
            kids.emplace_back(); taxv.push_back(t); solidv.push_back(1);
            int u = (int)kids.size() - 1;
            if (stk.empty()) root = u; else kids[stk.back()].push_back(u);
            i = j;
        } else { ++i; }
    }
    g.n = (int)kids.size();
    if (g.n == 0) return;
    g.head.assign(g.n + 1, 0);
    for (int u = 0; u < g.n; ++u) g.head[u + 1] = g.head[u] + (int)kids[u].size();
    g.ch.resize(g.head[g.n]);
    g.par.assign(g.n, -1);
    g.tax = taxv;
    g.solid.assign(g.n, 0);
    for (int u = 0; u < g.n; ++u) {
        if (!kids[u].empty())
            std::memcpy(&g.ch[g.head[u]], kids[u].data(), kids[u].size() * sizeof(int));
        g.solid[u] = (char)solidv[u];          // leaf edges are solid, by construction
        for (int c : kids[u]) g.par[c] = u;
    }
    g.nleaf = 0; g.maxdeg = 1;
    for (int u = 0; u < g.n; ++u) {
        g.nleaf += (g.tax[u] >= 0);
        g.maxdeg = std::max(g.maxdeg, g.head[u + 1] - g.head[u] + 1);
    }
    g.post.clear(); g.post.reserve(g.n);
    std::vector<int> st{root};
    while (!st.empty()) {
        int u = st.back(); st.pop_back();
        g.post.push_back(u);
        for (int k = g.head[u]; k < g.head[u + 1]; ++k) st.push_back(g.ch[k]);
    }
    std::reverse(g.post.begin(), g.post.end());
}

// species tree as a plain adjacency structure
struct AdjTree {
    int n = 0, nTaxa = 0;
    std::vector<int> head, nb, tax;
    std::vector<std::vector<int>> deg;
    int add(int t = -1) { deg.emplace_back(); tax.push_back(t); return n++; }
    void edge(int u, int v) { deg[u].push_back(v); deg[v].push_back(u); }
    void build() {
        head.assign(n + 1, 0);
        for (int u = 0; u < n; ++u) head[u + 1] = head[u] + (int)deg[u].size();
        nb.clear();
        for (int u = 0; u < n; ++u) nb.insert(nb.end(), deg[u].begin(), deg[u].end());
        nTaxa = 0;
        for (int u = 0; u < n; ++u) nTaxa = std::max(nTaxa, tax[u] + 1);
    }
};

static AdjTree parse_adj(const std::string &s)
{
    AdjTree t;
    std::vector<int> stk;
    size_t i = 0, n = s.size();
    while (i < n) {
        char c = s[i];
        if (c == '(') { int u = t.add(); stk.push_back(u); ++i; }
        else if (c == ',' || c == ';' || c == ' ' || c == '\t' || c == '\r' || c == '\n') ++i;
        else if (c == ')') {
            int u = stk.back(); stk.pop_back(); ++i;
            if (i < n && (s[i] == 'S' || s[i] == 'D')) ++i;
            if (!stk.empty()) t.edge(stk.back(), u);
        } else if (c >= '0' && c <= '9') {
            size_t j = i; int v = 0;
            while (j < n && s[j] >= '0' && s[j] <= '9') { v = v * 10 + (s[j] - '0'); ++j; }
            int u = t.add(v);
            if (!stk.empty()) t.edge(stk.back(), u);
            i = j;
        } else ++i;
    }
    bool changed = true;
    while (changed) {                       // suppress degree-2 nodes
        changed = false;
        for (int u = 0; u < t.n; ++u) {
            if (t.deg[u].size() == 2) {
                int a = t.deg[u][0], b = t.deg[u][1];
                t.deg[a].erase(std::find(t.deg[a].begin(), t.deg[a].end(), u));
                t.deg[b].erase(std::find(t.deg[b].begin(), t.deg[b].end(), u));
                t.edge(a, b);
                t.deg[u].clear();
                changed = true;
            }
        }
    }
    t.build();
    return t;
}

// --------------------------------------------------------------- utilities --
static std::string read_file(const char *path)
{
    FILE *f = std::fopen(path, "rb");
    if (!f) { std::fprintf(stderr, "cannot open %s\n", path); std::exit(1); }
    std::string s;
    char buf[1 << 20];
    size_t k;
    while ((k = std::fread(buf, 1, sizeof buf, f)) > 0) s.append(buf, k);
    std::fclose(f);
    return s;
}

static void split_lines(const std::string &s, std::vector<std::string> &out)
{
    size_t i = 0, n = s.size();
    while (i < n) {
        size_t j = s.find('\n', i);
        if (j == std::string::npos) j = n;
        if (j > i) out.emplace_back(s.substr(i, j - i));
        i = j + 1;
    }
}

// solid clusters: gamma = largest one (route A for its leaves),
// everything else gets route B (direct enumeration inside its own cluster).
static void mark_clusters(Gene &g)
{
    std::vector<int> uf(g.n);
    for (int u = 0; u < g.n; ++u) uf[u] = u;
    std::function<int(int)> find = [&](int a) {
        while (uf[a] != a) { uf[a] = uf[uf[a]]; a = uf[a]; }
        return a;
    };
    for (int u = 0; u < g.n; ++u)
        if (g.par[u] >= 0 && g.solid[u]) {
            int a = find(u), b = find(g.par[u]);
            if (a != b) uf[a] = b;
        }
    // index clusters, count leaves
    std::vector<int> idx(g.n, -1), cnt;
    for (int u = 0; u < g.n; ++u) {
        int r = find(u);
        if (idx[r] < 0) { idx[r] = (int)cnt.size(); cnt.push_back(0); }
    }
    g.ncl = (int)cnt.size();
    g.clId.assign(g.n, -1);
    for (int u = 0; u < g.n; ++u) {
        g.clId[u] = idx[find(u)];
        if (g.tax[u] >= 0) cnt[g.clId[u]]++;
    }
    int best = 0;
    for (int c = 0; c < g.ncl; ++c) if (cnt[c] > cnt[best]) best = c;
    g.bestCl = best;
    g.small.assign(g.n, 0);
    g.inGamma.assign(g.n, 0);
    for (int u = 0; u < g.n; ++u) {
        if (g.tax[u] < 0) continue;
        if (g.clId[u] == best) g.inGamma[u] = 1;
        else g.small[u] = 1;
    }
    // compact slot index: only nodes outside the dominant cluster need CPS storage
    g.slot.assign(g.n, -1);
    for (int u = 0; u < g.n; ++u)
        if (g.clId[u] != best) g.slot[u] = g.nslot++;
    // CSR of leaves per cluster
    g.clHead.assign(g.ncl + 1, 0);
    for (int u = 0; u < g.n; ++u)
        if (g.tax[u] >= 0) g.clHead[g.clId[u] + 1]++;
    for (int c = 0; c < g.ncl; ++c) g.clHead[c + 1] += g.clHead[c];
    g.clLeaf.resize(g.clHead[g.ncl]);
    std::vector<int> fill(g.clHead.begin(), g.clHead.end() - 1);
    for (int u = 0; u < g.n; ++u)
        if (g.tax[u] >= 0) g.clLeaf[fill[g.clId[u]]++] = u;
    // depth
    g.depth.assign(g.n, 0);
    for (int idx2 = g.n - 1; idx2 >= 0; --idx2) {
        int u = g.post[idx2];
        if (g.par[u] >= 0) g.depth[u] = g.depth[g.par[u]] + 1;
    }
}

// direction (1..3) of every taxon w.r.t. internal node m of the species tree
static void compute_sdir(const AdjTree &S, int m, uint8_t *sdir,
                         std::vector<int> &seen, int &stamp, std::vector<int> &stk)
{
    int k = (int)S.deg[m].size();
    for (int di = 0; di < k; ++di) {
        ++stamp;
        int start = S.deg[m][di];
        seen[m] = stamp;
        stk.clear(); stk.push_back(start);
        while (!stk.empty()) {
            int u = stk.back(); stk.pop_back();
            if (seen[u] == stamp) continue;
            seen[u] = stamp;
            if (S.tax[u] >= 0) sdir[S.tax[u]] = (uint8_t)(di + 1);
            for (int e = S.head[u]; e < S.head[u + 1]; ++e) {
                int w = S.nb[e];
                if (seen[w] != stamp) stk.push_back(w);
            }
        }
    }
}

// ------------------------------------------------------------------- sweep --
struct Scratch {
    std::vector<int> cntT, cntS, cntC, upT, upS, upC, Qv, Pv;
    std::vector<long long> prefP, prefC, E, Ec;
    std::vector<int> cpsUp, cpsDown;
    void ensure(int n, int maxdeg) {
        const size_t n3 = (size_t)3 * n;
        if (cntT.size() < n3) {
            cntT.resize(n3); cntS.resize(n3); cntC.resize(n3);
            upT.resize(n3); upS.resize(n3); upC.resize(n3);
            prefP.resize(n3); prefC.resize(n3);
            cpsUp.resize(n3); cpsDown.resize(n3);
        }
        const size_t k3 = (size_t)3 * maxdeg;
        if (Qv.size() < k3) { Qv.resize(k3); Pv.resize(k3); E.resize(k3); Ec.resize(k3); }
    }
};

// SUM over internal nodes q on path(u,v) of CPS_d(q, direction towards u)
static inline long long pathsum(const Gene &g, Scratch &sc, int u, int v, int d)
{
    const int n = g.n;
    const int *cUp = &sc.cpsUp[(size_t)d * n];
    const int *cDn = &sc.cpsDown[(size_t)d * n];
    long long s = 0;
    while (u != v) {
        if (g.depth[u] > g.depth[v]) { s += cDn[g.slot[u]]; u = g.par[u]; }
        else if (g.depth[v] > g.depth[u]) { if (g.tax[v] < 0) s += cUp[g.slot[v]]; v = g.par[v]; }
        else {
            s += cDn[g.slot[u]]; u = g.par[u];
            if (g.tax[v] < 0) s += cUp[g.slot[v]];   // a leaf contributes nothing
            v = g.par[v];
        }
    }
    return s;
}

static void sweep_gene(const Gene &g, const uint8_t *sdir, int nTaxa, Scratch &sc,
                       std::vector<long long> &accP, std::vector<long long> &accCPS)
{
    const int n = g.n;
    // cnt*/pref* are fully written below; only the root's up*/pref* need zeroing
    const int root0 = g.post.back();
    for (int d = 0; d < 3; ++d) {
        if (g.slot[root0] >= 0) {                 // root has no parent direction
            sc.cpsUp[(size_t)d * n + g.slot[root0]] = 0;
            sc.cpsDown[(size_t)d * n + g.slot[root0]] = 0;
        }
        sc.upT[(size_t)d * n + root0] = 0;
        sc.upS[(size_t)d * n + root0] = 0;
        sc.upC[(size_t)d * n + root0] = 0;
        sc.prefP[(size_t)d * n + root0] = 0;
        sc.prefC[(size_t)d * n + root0] = 0;
    }

    // ---- bottom-up
    for (int idx = 0; idx < n; ++idx) {
        const int u = g.post[idx];
        if (g.tax[u] >= 0) {
            const int j = sdir[g.tax[u]];
            for (int d = 0; d < 3; ++d) {
                const int val = (d + 1 == j) ? 1 : 0;
                sc.cntT[(size_t)d * n + u] = val;
                sc.cntS[(size_t)d * n + u] = val;
                sc.cntC[(size_t)d * n + u] = g.small[u] ? val : 0;
            }
        } else {
            for (int d = 0; d < 3; ++d) {
                int t = 0, s = 0, c = 0;
                for (int e = g.head[u]; e < g.head[u + 1]; ++e) {
                    const int w = g.ch[e];
                    t += sc.cntT[(size_t)d * n + w];
                    c += sc.cntC[(size_t)d * n + w];
                    if (g.solid[w]) s += sc.cntS[(size_t)d * n + w];
                }
                sc.cntT[(size_t)d * n + u] = t;
                sc.cntS[(size_t)d * n + u] = s;
                sc.cntC[(size_t)d * n + u] = c;
            }
        }
    }

    // ---- top-down
    long long totalP[3] = {0, 0, 0}, totalC[3] = {0, 0, 0};
    for (int idx = n - 1; idx >= 0; --idx) {
        const int u = g.post[idx];
        const int p = g.par[u];
        if (p >= 0) {
            for (int d = 0; d < 3; ++d) {
                const int allT = sc.upT[(size_t)d * n + p] + sc.cntT[(size_t)d * n + p];
                const int allS = sc.upS[(size_t)d * n + p] + sc.cntS[(size_t)d * n + p];
                const int allC = sc.upC[(size_t)d * n + p] + sc.cntC[(size_t)d * n + p];
                sc.upT[(size_t)d * n + u] = allT - sc.cntT[(size_t)d * n + u];
                sc.upC[(size_t)d * n + u] = allC - sc.cntC[(size_t)d * n + u];
                sc.upS[(size_t)d * n + u] = g.solid[u] ? (allS - sc.cntS[(size_t)d * n + u]) : 0;
            }
        }
        if (g.tax[u] >= 0) {
            if (g.inGamma[u]) {                     // route A
                const int t = g.tax[u];
                for (int d = 0; d < 3; ++d) {
                    accP[(size_t)d * nTaxa + t] += sc.prefP[(size_t)d * n + u];
                    accCPS[(size_t)d * nTaxa + t] += sc.prefC[(size_t)d * n + u];
                }
            }
            continue;                               // route B is done after the pass
        }
        const int nch = g.head[u + 1] - g.head[u];
        const int k = nch + (p >= 0 ? 1 : 0);

        for (int ci = 0; ci < nch; ++ci) {
            const int c = g.ch[g.head[u] + ci];
            for (int d = 0; d < 3; ++d) {
                sc.Qv[(size_t)d * k + ci] = g.solid[c] ? sc.cntS[(size_t)d * n + c] : 0;
                sc.Pv[(size_t)d * k + ci] = sc.cntT[(size_t)d * n + c] - sc.cntC[(size_t)d * n + c];
            }
        }
        if (p >= 0) {
            for (int d = 0; d < 3; ++d) {
                sc.Qv[(size_t)d * k + (k - 1)] = sc.upS[(size_t)d * n + u];
                sc.Pv[(size_t)d * k + (k - 1)] = sc.upT[(size_t)d * n + u] - sc.upC[(size_t)d * n + u];
            }
        }

        for (int d = 0; d < 3; ++d) {
            const int a = (d + 1) % 3, b = (d + 2) % 3;
            long long tA = 0, tB = 0, cross = 0;
            for (int i = 0; i < k; ++i) {
                const long long A = sc.Qv[(size_t)a * k + i], B = sc.Qv[(size_t)b * k + i];
                tA += A; tB += B; cross += A * B;
            }
            for (int i = 0; i < k; ++i) {
                const long long A = sc.Qv[(size_t)a * k + i], B = sc.Qv[(size_t)b * k + i];
                const long long cps = (tA - A) * (tB - B) - (cross - A * B);
                sc.E[(size_t)d * k + i] = sc.Pv[(size_t)d * k + i] * cps;
                sc.Ec[(size_t)d * k + i] = cps;
            }
            const long long gp = (p >= 0) ? sc.E[(size_t)d * k + (k - 1)] : 0;
            const long long gpc = (p >= 0) ? sc.Ec[(size_t)d * k + (k - 1)] : 0;
            for (int ci = 0; ci < nch; ++ci) {
                const int c = g.ch[g.head[u] + ci];
                if (g.slot[c] >= 0)
                    sc.cpsDown[(size_t)d * n + g.slot[c]] = (int)sc.Ec[(size_t)d * k + ci];
            }
            if (p >= 0) {
                totalP[d] += gp; totalC[d] += gpc;
                if (g.slot[u] >= 0)
                    sc.cpsUp[(size_t)d * n + g.slot[u]] = (int)gpc;
            }
            const long long baseP = sc.prefP[(size_t)d * n + u];
            const long long baseC = sc.prefC[(size_t)d * n + u];
            for (int ci = 0; ci < nch; ++ci) {
                const int c = g.ch[g.head[u] + ci];
                sc.prefP[(size_t)d * n + c] = baseP + sc.E[(size_t)d * k + ci] - gp;
                sc.prefC[(size_t)d * n + c] = baseC + sc.Ec[(size_t)d * k + ci] - gpc;
            }
        }
    }
    // the TOTALP constant of Theorem 3(10), route A leaves only
    for (int idx = 0; idx < n; ++idx) {
        const int u = g.post[idx];
        if (g.tax[u] < 0 || !g.inGamma[u]) continue;
        const int t = g.tax[u];
        for (int d = 0; d < 3; ++d) {
            accP[(size_t)d * nTaxa + t] += totalP[d];
            accCPS[(size_t)d * nTaxa + t] += totalC[d];
        }
    }
    // ---- route B: leaves outside the dominant cluster.
    // T_x(m,d) = SUM_{a in gamma_x ∩ D_d} ( S_x(d) - PathSum_x(a,d) )
    for (int idx = 0; idx < n; ++idx) {
        const int u = g.post[idx];
        if (g.tax[u] < 0 || g.inGamma[u]) continue;
        long long Sx[3], val[3] = {0, 0, 0};
        for (int d = 0; d < 3; ++d) Sx[d] = sc.prefC[(size_t)d * n + u] + totalC[d];
        const int cid = g.clId[u];
        for (int e = g.clHead[cid]; e < g.clHead[cid + 1]; ++e) {
            const int a = g.clLeaf[e];
            if (a == u) continue;
            const int d = sdir[g.tax[a]] - 1;
            if (d < 0 || d > 2) continue;
            val[d] += Sx[d] - pathsum(g, sc, u, a, d);
        }
        const int x = g.tax[u];
        for (int d = 0; d < 3; ++d) accP[(size_t)d * nTaxa + x] += val[d];
    }
}

// correction: T_x(m,c) -= S_x(m,c) for c = direction of x at m
static void apply_correction(std::vector<long long> &accP, const std::vector<long long> &accCPS,
                             const uint8_t *sdir, int nTaxa)
{
    for (int t = 0; t < nTaxa; ++t) {
        const int c = sdir[t];
        if (c < 1 || c > 3) continue;
        accP[(size_t)(c - 1) * nTaxa + t] -= accCPS[(size_t)(c - 1) * nTaxa + t];
    }
}

int main(int argc, char **argv)
{
    const char *spPath = nullptr, *gnPath = nullptr;
    std::string strategy = "m";
    int threads = omp_get_max_threads();
    int mstart = 0, mcount = 1000000, block = 8;
    bool doStore = false;
    const char *dumpPath = nullptr;
    const char *dumpBinPath = nullptr;

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto nxt = [&]() { return std::string(argv[++i]); };
        if (a == "--species") spPath = strdup(nxt().c_str());
        else if (a == "--genes") gnPath = strdup(nxt().c_str());
        else if (a == "--strategy") strategy = nxt();
        else if (a == "--threads") threads = atoi(nxt().c_str());
        else if (a == "--block") block = atoi(nxt().c_str());
        else if (a == "--m-start") mstart = atoi(nxt().c_str());
        else if (a == "--m-count") mcount = atoi(nxt().c_str());
        else if (a == "--store") doStore = true;
        else if (a == "--dump") dumpPath = strdup(nxt().c_str());
        else if (a == "--dumpbin") dumpBinPath = strdup(nxt().c_str());
        else { std::fprintf(stderr, "unknown arg %s\n", a.c_str()); return 1; }
    }
    if (!spPath || !gnPath) { std::fprintf(stderr, "need --species and --genes\n"); return 1; }

    AdjTree S = parse_adj(read_file(spPath));
    std::vector<std::string> lines;
    split_lines(read_file(gnPath), lines);
    std::vector<Gene> genes;
    genes.reserve(lines.size());
    size_t K = 0;
    for (auto &ln : lines) {
        if (ln.find('(') == std::string::npos) continue;
        Gene g;
        parse_newick(ln, g);
        if (g.n < 3) continue;
        mark_clusters(g);
        K += g.n;
        genes.push_back(std::move(g));
    }
    const int nTaxa = S.nTaxa;
    std::vector<int> ms;
    for (int u = 0; u < S.n && (int)ms.size() < mstart + mcount; ++u)
        if (S.deg[u].size() >= 3) ms.push_back(u);
    if (mstart > 0) ms.erase(ms.begin(), ms.begin() + std::min<size_t>(mstart, ms.size()));
    if ((int)ms.size() > mcount) ms.resize(mcount);
    std::fprintf(stderr, "[sector] species: n=%d taxa=%d | genes=%zu K=%zu | m's=%zu | threads=%d | strategy=%s",
                 S.n, nTaxa, genes.size(), K, ms.size(), threads, strategy.c_str());
    if (strategy == "block") std::fprintf(stderr, " (B=%d)", block);
    std::fprintf(stderr, "\n");

    std::vector<long long> store;
    if (doStore) store.assign((size_t)3 * nTaxa * ms.size(), 0);

    FILE *df = nullptr;
    if (dumpPath) {
        df = std::fopen(dumpPath, "w");
        if (!df) { std::fprintf(stderr, "cannot write %s\n", dumpPath); return 1; }
        if (threads != 1) { std::fprintf(stderr, "[sector] --dump forces 1 thread\n"); threads = 1; }
    }
    FILE *bf = nullptr;
    if (dumpBinPath) {
        bf = std::fopen(dumpBinPath, "wb");
        if (!bf) { std::fprintf(stderr, "cannot write %s\n", dumpBinPath); return 1; }
        // header: M, nTaxa
        const long long hdr[2] = {(long long)ms.size(), (long long)nTaxa};
        std::fwrite(hdr, sizeof(long long), 2, bf);
        std::fflush(bf);
    }
    omp_set_num_threads(threads);
    const double t0 = omp_get_wtime();
    long long checksum = 0, negcount = 0;

    if (strategy == "m" || strategy == "block") {
        const int B = (strategy == "block") ? block : 1;
        const int nblk = (int)((ms.size() + B - 1) / B);
#pragma omp parallel
        {
            Scratch sc;
            std::vector<int> seen(S.n, -1), stk;
            stk.reserve(S.n);
            int stamp = 0;
            std::vector<uint8_t> dirs((size_t)B * nTaxa);
            std::vector<std::vector<long long>> accP(B), accCPS(B);
#pragma omp for schedule(dynamic, 1)
            for (int bi = 0; bi < nblk; ++bi) {
                const int lo = bi * B, hi = std::min<int>((int)ms.size(), lo + B);
                const int nb = hi - lo;
                for (int b = 0; b < nb; ++b) {
                    compute_sdir(S, ms[lo + b], &dirs[(size_t)b * nTaxa], seen, stamp, stk);
                    accP[b].assign((size_t)3 * nTaxa, 0LL);
                    accCPS[b].assign((size_t)3 * nTaxa, 0LL);
                }
                // gene-major: the gene structure stays hot across the whole block
                for (auto &g : genes) { sc.ensure(g.n, g.maxdeg);
                    for (int b = 0; b < nb; ++b)
                        sweep_gene(g, &dirs[(size_t)b * nTaxa], nTaxa, sc, accP[b], accCPS[b]);
                }
                long long s = 0;
                for (int b = 0; b < nb; ++b) {
                    apply_correction(accP[b], accCPS[b], &dirs[(size_t)b * nTaxa], nTaxa);
                    if (bf) {
                        const long long off = (long long)(16 + ((size_t)(lo + b) * 3 * nTaxa) * sizeof(long long));
#pragma omp critical(writebin)
                        {
                            std::fseek(bf, (long)off, SEEK_SET);
                            std::fwrite(accP[b].data(), sizeof(long long), (size_t)3 * nTaxa, bf);
                        }
                    }
                    if (dumpPath) {
                        long long c1 = 0, c2 = 0, c3 = 0;
                        for (int t = 0; t < nTaxa; ++t) {
                            const int q = dirs[(size_t)b * nTaxa + t];
                            c1 += (q == 1); c2 += (q == 2); c3 += (q == 3);
                        }
                        std::fprintf(df, "M %d %lld %lld %lld\n", ms[lo + b], c1, c2, c3);
                        for (int t = 0; t < nTaxa; ++t)
                            std::fprintf(df, "T %d %d %lld %lld %lld\n", ms[lo + b], t,
                                         accP[b][(size_t)0 * nTaxa + t],
                                         accP[b][(size_t)1 * nTaxa + t],
                                         accP[b][(size_t)2 * nTaxa + t]);
                    }
                    if (doStore)
                        std::memcpy(&store[(size_t)(lo + b) * 3 * nTaxa], accP[b].data(),
                                    (size_t)3 * nTaxa * sizeof(long long));
                    for (auto v : accP[b]) { s += v; if (v < 0) ++negcount; }
                }
#pragma omp atomic
                checksum += s;
            }
        }
    } else if (strategy == "nested") {
        // all threads cooperate on one m: per-thread accumulators + reduction per m
        const int G = (int)genes.size();
        // these must be SHARED across the team
        std::vector<std::vector<long long>> allP(threads), allC(threads);
        for (int t = 0; t < threads; ++t) {
            allP[t].assign((size_t)3 * nTaxa, 0LL);
            allC[t].assign((size_t)3 * nTaxa, 0LL);
        }
        std::vector<long long> redP((size_t)3 * nTaxa), redC((size_t)3 * nTaxa);
        std::vector<uint8_t> sd(nTaxa);
        int stamp = 0;
#pragma omp parallel
        {
            const int tid = omp_get_thread_num(), nth = omp_get_num_threads();
            Scratch sc;
            std::vector<int> seen(S.n, -1), stk;
            stk.reserve(S.n);
            for (size_t mi = 0; mi < ms.size(); ++mi) {
#pragma omp single
                compute_sdir(S, ms[mi], sd.data(), seen, stamp, stk);
                std::fill(allP[tid].begin(), allP[tid].end(), 0LL);
                std::fill(allC[tid].begin(), allC[tid].end(), 0LL);
                for (int gi = tid; gi < G; gi += nth) {
                    sc.ensure(genes[gi].n, genes[gi].maxdeg);
                    sweep_gene(genes[gi], sd.data(), nTaxa, sc, allP[tid], allC[tid]);
                }
#pragma omp barrier
#pragma omp single
                {
                    std::fill(redP.begin(), redP.end(), 0LL);
                    std::fill(redC.begin(), redC.end(), 0LL);
                    for (int t = 0; t < nth; ++t)
                        for (size_t z = 0; z < redP.size(); ++z) {
                            redP[z] += allP[t][z];
                            redC[z] += allC[t][z];
                        }
                    apply_correction(redP, redC, sd.data(), nTaxa);
                    long long s = 0;
                    for (auto v : redP) { s += v; if (v < 0) ++negcount; }
                    checksum += s;
                }
            }
        }
    } else {
        std::fprintf(stderr, "unknown strategy\n");
        return 1;
    }

    if (df) std::fclose(df);
    if (bf) std::fclose(bf);
    const double dt = omp_get_wtime() - t0;
    const double perm = (double)ms.size();
    std::printf("RESULT strategy=%s threads=%d m=%zu block=%d time=%.3f s  %.1f ms/m  "
                "%.3e nodevisits/s  checksum=%lld neg=%lld\n",
                strategy.c_str(), threads, ms.size(), block, dt, 1000.0 * dt / perm,
                (double)K * perm / dt, checksum, negcount);
    return 0;
}
