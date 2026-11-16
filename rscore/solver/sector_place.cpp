// sector_place.cpp — phase 2 of the sector method.
//
// Reads the T table produced by
//     sector_sweep --strategy m --dump T.txt ...
// and, for every leaf x, computes the exact score change of re-attaching x at every
// edge of S-x (Theorem 2: rerooting, only the crossed node's direction changes):
//
//     G_x(R) = SUM over internal m of S, m != parent(x) :  T_x(m, dir_m(R))
//     Delta(x, R) = G_x(R) - G_x(R0),   R0 = the position x currently occupies
//
// It then prints the ranked moves and optionally applies the best N of them.
//
// Build: g++ -std=c++20 -O3 -march=native -fopenmp sector_place.cpp -o sector_place

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <functional>
#include <string>
#include <vector>

struct AdjTree {
    int n = 0, nTaxa = 0;
    std::vector<int> head, nb, tax, tail;
    std::vector<std::vector<int>> deg;
    int add(int t = -1) { deg.emplace_back(); tax.push_back(t); return n++; }
    void edge(int u, int v) { deg[u].push_back(v); deg[v].push_back(u); }
    void build() {
        head.assign(n + 1, 0);
        for (int u = 0; u < n; ++u) head[u + 1] = head[u] + (int)deg[u].size();
        nb.clear(); tail.clear();
        for (int u = 0; u < n; ++u) {
            nb.insert(nb.end(), deg[u].begin(), deg[u].end());
            tail.insert(tail.end(), deg[u].size(), u);
        }
        nTaxa = 0;
        for (int u = 0; u < n; ++u) nTaxa = std::max(nTaxa, tax[u] + 1);
    }
    int pos(int u, int v) const {                  // CSR position of directed edge u->v
        for (int e = head[u]; e < head[u + 1]; ++e) if (nb[e] == v) return e;
        return -1;
    }
    int degOf(int u) const { return head[u + 1] - head[u]; }
};

static std::string read_file(const char *path)
{
    FILE *f = std::fopen(path, "rb");
    if (!f) { std::fprintf(stderr, "cannot open %s\n", path); std::exit(1); }
    std::string s; char buf[1 << 20]; size_t k;
    while ((k = std::fread(buf, 1, sizeof buf, f)) > 0) s.append(buf, k);
    std::fclose(f);
    return s;
}

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
    while (changed) {
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

static void write_newick(const AdjTree &t, int u, int p, std::string &out)
{
    if (t.tax[u] >= 0) { out += std::to_string(t.tax[u]); return; }
    out += '(';
    bool first = true;
    for (int v : t.deg[u]) {
        if (v == p) continue;
        if (!first) out += ',';
        first = false;
        write_newick(t, v, u, out);
    }
    out += ')';
}

int main(int argc, char **argv)
{
    const char *spPath = nullptr, *tPath = nullptr, *outPath = nullptr, *tbinPath = nullptr;
    int applyN = 0, onlyLeaf = -1;
    bool dumpAll = false;
    std::string probeSpec;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto nxt = [&]() { return std::string(argv[++i]); };
        if (a == "--species") spPath = strdup(nxt().c_str());
        else if (a == "--T") tPath = strdup(nxt().c_str());
        else if (a == "--Tbin") tbinPath = strdup(nxt().c_str());
        else if (a == "--out") outPath = strdup(nxt().c_str());
        else if (a == "--apply") applyN = atoi(nxt().c_str());
        else if (a == "--leaf") onlyLeaf = atoi(nxt().c_str());
        else if (a == "--all") dumpAll = true;
        else if (a == "--probe") probeSpec = nxt();
        else { std::fprintf(stderr, "unknown arg %s\n", a.c_str()); return 1; }
    }
    if (!spPath || (!tPath && !tbinPath)) { std::fprintf(stderr, "need --species and --T/--Tbin\n"); return 1; }

    AdjTree S = parse_adj(read_file(spPath));
    const int nTaxa = S.nTaxa;
    int nTaxa2 = 0;
    for (int q = 0; q < S.n; ++q) nTaxa2 += (S.tax[q] >= 0);

    // internal nodes with >= 3 directions, ascending -- the order sector_sweep used
    std::vector<int> ms;
    for (int u = 0; u < S.n; ++u) if (S.deg[u].size() >= 3) ms.push_back(u);
    const int M = (int)ms.size();

    // ---- read the dump into T[mi][taxon][d]
    std::vector<long long> T((size_t)M * nTaxa * 3, 0);
    if (tbinPath) {
        FILE *f = std::fopen(tbinPath, "rb");
        if (!f) { std::fprintf(stderr, "cannot open %s\n", tbinPath); return 1; }
        long long hdr[2] = {0, 0};
        if (std::fread(hdr, sizeof(long long), 2, f) != 2) { std::fprintf(stderr, "bad header\n"); return 1; }
        if (hdr[0] != M || hdr[1] != nTaxa) {
            std::fprintf(stderr, "T size mismatch: file M=%lld nTaxa=%lld, here M=%d nTaxa=%d\n",
                         hdr[0], hdr[1], M, nTaxa);
            return 1;
        }
        // file layout is [m][d][taxon]; T is [m][taxon][d] -> transpose
        std::vector<long long> blk((size_t)3 * nTaxa);
        for (int mi = 0; mi < M; ++mi) {
            if (std::fread(blk.data(), sizeof(long long), blk.size(), f) != blk.size()) {
                std::fprintf(stderr, "short read at m-index %d\n", mi); return 1;
            }
            for (int d = 0; d < 3; ++d)
                for (int t = 0; t < nTaxa; ++t)
                    T[((size_t)mi * nTaxa + t) * 3 + d] = blk[(size_t)d * nTaxa + t];
        }
        std::fclose(f);
        std::fprintf(stderr, "[place] read %zu T entries from %s\n", T.size(), tbinPath);
    }
    std::vector<int> midx(S.n, -1);
    for (int i = 0; i < M; ++i) midx[ms[i]] = i;
    if (!tbinPath) {
        std::string str = read_file(tPath);
        const char *b = str.c_str(), *e = b + str.size();
        const char *q = b;
        while (q < e) {
            const char *nl = (const char *)std::memchr(q, '\n', (size_t)(e - q));
            if (!nl) nl = e;
            if (nl > q + 4 && *q == 'T') {
                const char *r = q + 1;
                char *endp = nullptr;
                long m = std::strtol(r, &endp, 10); r = endp;
                long tax = std::strtol(r, &endp, 10); r = endp;
                long long v0 = std::strtoll(r, &endp, 10); r = endp;
                long long v1 = std::strtoll(r, &endp, 10); r = endp;
                long long v2 = std::strtoll(r, &endp, 10);
                if (m >= 0 && m < S.n && tax >= 0 && tax < nTaxa) {
                    int mi = midx[m];
                    if (mi >= 0) {
                        size_t base = ((size_t)mi * nTaxa + tax) * 3;
                        T[base + 0] = v0; T[base + 1] = v1; T[base + 2] = v2;
                    }
                }
            }
            q = nl + 1;
        }
    }
    std::fprintf(stderr, "[place] T loaded, M=%d nTaxa=%d\n", M, nTaxa);
    // ---- phase 2: exact Delta for every (leaf, edge of S-x)
    struct Move { long long delta; int x, u, v; };
    std::vector<Move> moves;
    std::vector<long long> G((size_t)S.nb.size(), 0);
    std::vector<int> bfsOrder;

    for (int xi = 0; xi < S.n; ++xi) {
        const int x = xi;
        if (S.tax[x] < 0) continue;
        if (S.deg[x].size() != 1) continue;
        if (onlyLeaf >= 0 && S.tax[x] != onlyLeaf) continue;
        const int p = S.deg[x][0];
        if (S.deg[p].size() != 3) continue;            // species tree must be binary
        int ab[2], nab = 0;
        for (int v : S.deg[p]) if (v != x) ab[nab++] = v;
        if (nab != 2) continue;

        // dirOf(m, target) = index of the direction of m whose component contains `target`.
        // (Was a bug: checking adjacency instead of component membership.)
        auto hopTable = [&](int root, std::vector<int> &par) {
            // BFS from `root`: par[z] = the neighbour of z on the path towards root
            par.assign(S.n, -1);
            std::vector<int> q{root};
            par[root] = root;
            for (size_t i = 0; i < q.size(); ++i) {
                const int z = q[i];
                for (int e = S.head[z]; e < S.head[z + 1]; ++e) {
                    const int w = S.nb[e];
                    if (par[w] < 0) { par[w] = z; q.push_back(w); }
                }
            }
        };
        auto dirIdxOfNbr = [&](int m, int w) {
            for (int e = S.head[m]; e < S.head[m + 1]; ++e)
                if (S.nb[e] == w) return e - S.head[m];
            return -1;
        };
        // G at x's current position: for every m, direction = component containing p/a
        std::vector<int> parP;
        hopTable(p, parP);
        long long G0 = 0;
        for (int m : ms) {
            if (m == p) continue;
            const int hop = parP[m];                    // neighbour of m towards p
            const int d = dirIdxOfNbr(m, hop);
            if (d >= 0) G0 += T[((size_t)midx[m] * nTaxa + S.tax[x]) * 3 + d];
        }
        // rerooting: moving the root point across a node changes only that node's term
        auto tval = [&](int m, int dir) -> long long {
            if (m == p || midx[m] < 0 || dir < 0) return 0;   // p is not a valid m
            return T[((size_t)midx[m] * nTaxa + S.tax[x]) * 3 + dir];
        };
        std::vector<long long> G(S.nb.size(), LLONG_MIN);
        std::vector<int> bfsOrder;
        const int startPos = S.pos(p, ab[0]);
        if (startPos >= 0) {
            const int rev = S.pos(ab[0], p);
            G[startPos] = G0;
            bfsOrder.push_back(startPos);
            if (rev >= 0) { G[rev] = G0; bfsOrder.push_back(rev); }
            std::vector<char> seen(S.nb.size(), 0);
            seen[startPos] = 1;
            if (rev >= 0) seen[rev] = 1;
            for (size_t qi = 0; qi < bfsOrder.size(); ++qi) {
                const int e = bfsOrder[qi];
                const int u = S.tail[e], v = S.nb[e];
                const long long g = G[e];
                if (midx[v] < 0) continue;                 // a leaf cannot be crossed
                const long long tOld = tval(v, dirIdxOfNbr(v, u));
                for (int f = S.head[v]; f < S.head[v + 1]; ++f) {
                    const int w = S.nb[f];
                    if (w == u || seen[f]) continue;
                    G[f] = g - tOld + tval(v, dirIdxOfNbr(v, w));
                    if (midx[w] >= 0) { seen[f] = 1; bfsOrder.push_back(f); }
                }
            }
        }
        long long best = (onlyLeaf >= 0 || dumpAll) ? (long long)-4e18 : 0;
        int bu = -1, bv = -1;
        for (int u = 0; u < S.n; ++u) {
            if (u == x || u == p) continue;
            for (int e = S.head[u]; e < S.head[u + 1]; ++e) {
                const int v = S.nb[e];
                if (v == x || v == p || v < u) continue;
                long long gv = std::max(G[S.pos(u, v)], G[S.pos(v, u)]);
                if (gv == LLONG_MIN) continue;
                const long long d = gv - G0;
                if (d > best) { best = d; bu = u; bv = v; }
            }
        }
        if (bu >= 0 && dumpAll) std::printf("A %d %lld\n", S.tax[x], best);
        if (bu >= 0) { moves.push_back({best, x, bu, bv});
            if (onlyLeaf >= 0) std::printf("[probe] leaf %d best edge (%d,%d) delta=%lld\n", S.tax[x], bu, bv, best); }
        if ((xi % 100) == 0) std::fprintf(stderr, "[place] leaf %d/%d\n", xi, S.n);
    }

    std::sort(moves.begin(), moves.end(), [](const Move &a, const Move &b) { return a.delta > b.delta; });
    std::printf("[place] %zu leaves with a strictly improving move; best 10:\n", moves.size());
    for (int i = 0; i < (int)moves.size() && i < 10; ++i)
        std::printf("   +%lld  leaf %d -> edge (%d,%d)\n", moves[i].delta, S.tax[moves[i].x],
                    S.tax[moves[i].u] >= 0 ? S.tax[moves[i].u] : moves[i].u,
                    S.tax[moves[i].v] >= 0 ? S.tax[moves[i].v] : moves[i].v);
    long long predsum = 0;
    for (int i = 0; i < (int)moves.size() && i < applyN; ++i) predsum += moves[i].delta;
    std::printf("[place] applying %d moves, predicted gain = %lld\n", std::min<int>(applyN, (int)moves.size()), predsum);

    // ---- apply the chosen moves
    int applied = 0;
    std::vector<char> used(S.n, 0);
    for (auto &mv : moves) {
        if (applied >= applyN) break;
        const int x = mv.x, u = mv.u, v = mv.v;
        if (used[x] || used[u] || used[v]) continue;          // one move per leaf / edge
        const int p = S.deg[x][0];
        if (S.deg[x].size() != 1 || S.deg[p].size() != 3) continue;
        if (S.pos(u, v) < 0) continue;
        AdjTree const backup = S;                       // revert if the surgery breaks binarity
        int a = -1, b = -1;
        for (int z : S.deg[p]) if (z != x) { if (a < 0) a = z; else b = z; }
        if (a < 0 || b < 0) continue;
        auto drop = [](std::vector<int> &w, int z) {
            w.erase(std::remove(w.begin(), w.end(), z), w.end());
        };
        // detach x: suppress p, joining a--b
        drop(S.deg[a], p); drop(S.deg[b], p);
        drop(S.deg[p], a); drop(S.deg[p], b); drop(S.deg[p], x);
        S.deg[x].clear();
        S.edge(a, b);
        // re-attach x in the middle of edge (u,v)
        drop(S.deg[u], v); drop(S.deg[v], u);
        const int z = S.add();
        S.edge(z, u); S.edge(z, v); S.edge(z, x);
        bool ok = true;
        for (int q = 0; q < S.n; ++q) if (S.deg[q].size() > 3) { ok = false; break; }
        if (ok) {                                       // every internal node must stay <=3
            int nl = 0;
            for (int q = 0; q < S.n; ++q) nl += (S.tax[q] >= 0);
            if (nl != nTaxa2) ok = false;
        }
        if (!ok) { S = backup; continue; }
        used[x] = 1;
        ++applied;
        S.build();
    }
    std::printf("[place] applied %d moves\n", applied);

    if (outPath) {
        std::string out;
        int root = -1;
        for (int u = 0; u < S.n; ++u) if (S.tax[u] < 0 && S.deg[u].size() == 3) { root = u; break; }
        // root on an EDGE so the newick has a bifurcating root (scorer requires binary)
        const int a = S.deg[root][0], b = S.deg[root][1], c = S.deg[root][2];
        out += "(";
        write_newick(S, a, root, out);
        out += ",(";
        write_newick(S, b, root, out);
        out += ",";
        write_newick(S, c, root, out);
        out += "))";
        out += ";\n";
        FILE *f = std::fopen(outPath, "w");
        std::fputs(out.c_str(), f);
        std::fclose(f);
        std::printf("[place] wrote %s\n", outPath);
    }
    return 0;
}
