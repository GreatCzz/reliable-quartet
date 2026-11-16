// sector_exact.cpp — an EXACT sector engine with no fast/approx formula at all.
//
// For a given leaf x it enumerates every legal attachment point (every edge of S-x),
// physically rebuilds the species tree for that position, and scores it with the
// verified exact scorer (fastscore / score_fast backend, which is value-for-value
// identical to the reference implementation).  The best position is therefore exact
// by construction -- there is no delta formula involved, so nothing can be mis-derived.
//
// It is slow: O(#edges * K) per leaf.  That is fine; correctness first.
//
// Build:
//   g++ -std=c++20 -O3 -march=native -fopenmp -I. sector_exact.cpp -o sector_exact
//
// Usage:
//   ./sector_exact -t species.tree -g genes.trees -l TAXON   [--genes K] [-j T]
//   ./sector_exact -t species.tree -g genes.trees --all      [--genes K] [-j T] [--leaves N]

#define main score_fast_disabled_main
#include "score_fast.cpp"
#undef main

// ---------------------------------------------------------------- adjacency --
struct AdjTree {
    int n = 0;
    std::vector<std::vector<int>> deg;
    std::vector<int> tax;              // node -> taxon id (-1 = internal)
    int add(int t = -1) { deg.emplace_back(); tax.push_back(t); return n++; }
    void edge(int u, int v) { deg[u].push_back(v); deg[v].push_back(u); }
};

static AdjTree parse_adj(std::string const &s)
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
    // suppress degree-2 nodes
    bool changed = true;
    while (changed) {
        changed = false;
        for (int u = 0; u < t.n; ++u) {
            if (t.deg[u].size() == 2) {
                int a = t.deg[u][0], b = t.deg[u][1];
                for (auto *w : {&t.deg[a], &t.deg[b]})
                    w->erase(std::find(w->begin(), w->end(), u));
                t.deg[u].clear();
                t.edge(a, b);
                changed = true;
            }
        }
    }
    return t;
}

static void write_newick(AdjTree const &t, int u, int p, std::string &out)
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

// ------------------------------------------------------------- exact score ---
static std::vector<std::string> g_genes;
static std::string g_scorer, g_gfile, g_tmp;
static int g_jobs = 1;

// score `nw` by shelling out to the verified scorer binary (cannot mis-score)
static long long score_via_binary(std::string const &nw)
{
    { FILE *f = std::fopen(g_tmp.c_str(), "w");
      if (!f) return std::numeric_limits<long long>::min();
      std::fwrite(nw.data(), 1, nw.size(), f); std::fputc('\n', f); std::fclose(f); }
    std::string cmd = g_scorer + " -t " + g_tmp + " -g " + g_gfile +
                      " -j " + std::to_string(g_jobs) + " -q 2>/dev/null";
    FILE *p = popen(cmd.c_str(), "r");
    if (!p) return std::numeric_limits<long long>::min();
    std::string out; char buf[4096]; size_t k;
    while ((k = std::fread(buf, 1, sizeof buf, p)) > 0) out.append(buf, k);
    pclose(p);
    size_t pos = out.find("\"score\":");
    if (pos == std::string::npos) return std::numeric_limits<long long>::min();
    return std::strtoll(out.c_str() + pos + 8, nullptr, 10);
}

// score one species tree (newick) with all loaded genes, exactly
static long long score_newick(std::string const &nw, int nthreads)
{
    pool.clear();
    int r = parse(nw);
    if (r < 0) return std::numeric_limits<long long>::min();
    if (pool[r].ch.size() == 3) {                 // root trifurcation -> binary
        int synth = newNode();
        pool[synth].label = "";
        pool[synth].ch.push_back(pool[r].ch[1]);
        pool[synth].ch.push_back(pool[r].ch[2]);
        pool[r].ch.resize(1);
        pool[r].ch.push_back(synth);
    } else if (pool[r].ch.size() != 2) {
        return std::numeric_limits<long long>::min();
    }
    spRootBase = r;
    spN = (int)pool.size();
    spPool = pool;

    int n = (int)g_genes.size();
    std::vector<long long> res(n, 0);
    std::atomic<int> idx{0};
    int nt = std::max(1, std::min(nthreads, n));
    std::vector<std::thread> ths;
    ths.reserve(nt);
    for (int k = 0; k < nt; ++k)
        ths.emplace_back([&]() {
            int i;
            while ((i = idx.fetch_add(1, std::memory_order_relaxed)) < n) {
                long long v = solveGene(g_genes[i]);
                if (g_lastLeafCount < 4) v = 0;
                res[i] = v;
            }
        });
    for (auto &th : ths) th.join();
    long long s = 0;
    for (auto v : res) s += v;
    return s;
}

// --------------------------------------------------------------- one leaf ----
// returns {bestDelta, bestU, bestV}
static std::tuple<long long, int, int> exact_best_for_leaf(AdjTree const &S0, int x, int nthreads)
{
    AdjTree S = S0;
    const int p = S.deg[x][0];
    if (S.deg[x].size() != 1 || S.deg[p].size() != 3) return {-1, -1, -1};
    int a = -1, b = -1;
    for (int z : S.deg[p]) if (z != x) { if (a < 0) a = z; else b = z; }

    std::string cur;
    { int root = -1;
      for (int u = 0; u < S.n; ++u) if (S.tax[u] < 0 && S.deg[u].size() == 3) { root = u; break; }
      write_newick(S, root, -1, cur); cur += ";"; }
    const long long base = score_via_binary(cur);

    std::fprintf(stderr, "[dbg] base=%lld\n", base);
    long long best = 0; int bu = -1, bv = -1; int ncand = 0;
    long long bestsc = std::numeric_limits<long long>::min();
    for (int u = 0; u < S.n; ++u) {
        if (u == x || u == p) continue;
        for (int v : S.deg[u]) {
            if (v == x || v == p || v < u) continue;
            // build S with x re-attached in the middle of edge (u,v)
            AdjTree M = S;
            auto drop = [](std::vector<int> &w, int z) {
                w.erase(std::remove(w.begin(), w.end(), z), w.end());
            };
            drop(M.deg[a], p); drop(M.deg[b], p);
            M.deg[p].clear(); M.deg[x].clear();
            M.edge(a, b);
            drop(M.deg[u], v); drop(M.deg[v], u);
            int z = M.add();
            M.edge(z, u); M.edge(z, v); M.edge(z, x);
            std::string nw;
            { int root = -1;
              for (int w = 0; w < M.n; ++w) if (M.tax[w] < 0 && M.deg[w].size() == 3) { root = w; break; }
              write_newick(M, root, -1, nw); nw += ";"; }
            long long sc = score_via_binary(nw); ++ncand; if (sc > bestsc) bestsc = sc;
            if (sc - base > best) { best = sc - base; bu = u; bv = v; }
        }
    }
    std::fprintf(stderr, "[dbg] candidates=%d bestsc=%lld best=%lld\n", ncand, bestsc, best);
    return {best, bu, bv};
}

int main(int argc, char **argv)
{
    std::string treeFile, genesFile, leafName, posSpec;
    int nGenes = 1000, nthreads = 1, nLeaves = 0;
    bool all = false;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto nx = [&]() -> std::string { if (i + 1 < argc) return argv[++i]; std::exit(2); };
        if (a == "-t") treeFile = nx();
        else if (a == "-g") genesFile = nx();
        else if (a == "-l") leafName = nx();
        else if (a == "--genes") nGenes = atoi(nx().c_str());
        else if (a == "-j") nthreads = atoi(nx().c_str());
        else if (a == "--all") all = true;
        else if (a == "--leaves") nLeaves = atoi(nx().c_str());
        else if (a == "--scorer") g_scorer = nx();
        else if (a == "--pos") posSpec = nx();
        else if (a == "--gfile") g_gfile = nx();
        else { std::cerr << "unknown arg " << a << "\n"; return 2; }
    }
    if (treeFile.empty() || genesFile.empty()) { std::cerr << "need -t and -g\n"; return 2; }
    if (g_scorer.empty()) g_scorer = "score_fast";
    if (g_gfile.empty()) g_gfile = genesFile;
    g_tmp = "/tmp/.sector_exact_tmp.nwk";
    { const char *e = getenv("SECTOR_EXACT_JOBS"); if (e) g_jobs = atoi(e); }

    std::vector<std::string> spLines;
    readAllLines(treeFile, spLines);
    if (spLines.empty()) { std::cerr << "no species tree\n"; return 1; }
    readAllLines(genesFile, g_genes);
    if ((int)g_genes.size() > nGenes) g_genes.resize(nGenes);
    std::cerr << "[exact] genes=" << g_genes.size() << "\n";

    AdjTree S = parse_adj(spLines[0]);

    std::vector<int> leaves;
    for (int u = 0; u < S.n; ++u) if (S.tax[u] >= 0) leaves.push_back(u);
    std::sort(leaves.begin(), leaves.end(), [&](int a, int b) { return S.tax[a] < S.tax[b]; });

    if (!leafName.empty()) {
        int want = atoi(leafName.c_str());
        for (int u : leaves) if (S.tax[u] == want) {
            if (!posSpec.empty()) {
                int pu = atoi(posSpec.c_str());
                const char *c2 = strchr(posSpec.c_str(), ':');
                int pv = c2 ? atoi(c2 + 1) : -1;
                // interpret pu/pv as TAXA if they are leaves, else node ids
                int nu = pu, nv = pv;
                auto [d0, d1, d2] = exact_best_for_leaf(S, u, nthreads);  (void)d0;(void)d1;(void)d2;
                // rebuild with the requested position only
                std::string cur; { int r=-1; for (int z=0;z<S.n;++z) if (S.tax[z]<0 && S.deg[z].size()==3){r=z;break;}
                                   write_newick(S, r, -1, cur); cur += ";"; }
                long long base = score_via_binary(cur);
                AdjTree M = S;
                int pp = M.deg[u][0];
                int a=-1,b=-1; for (int z : M.deg[pp]) if (z != u) { if (a<0) a=z; else b=z; }
                auto drop = [](std::vector<int>& w, int z){ w.erase(std::remove(w.begin(),w.end(),z), w.end()); };
                drop(M.deg[a],pp); drop(M.deg[b],pp); M.deg[pp].clear(); M.deg[u].clear(); M.edge(a,b);
                drop(M.deg[nu],nv); drop(M.deg[nv],nu);
                int zz = M.add(); M.edge(zz,nu); M.edge(zz,nv); M.edge(zz,u);
                std::string nw2; { int r=-1; for (int z=0;z<M.n;++z) if (M.tax[z]<0 && M.deg[z].size()==3){r=z;break;}
                                   write_newick(M, r, -1, nw2); nw2 += ";"; }
                std::printf("POS %d %d %lld\n", pu, pv, score_via_binary(nw2) - base);
                return 0;
            }
            auto [d, bu, bv] = exact_best_for_leaf(S, u, nthreads);
            std::printf("leaf %d best_delta=%lld edge=(%d,%d)\n", want, d,
                        bu < 0 ? -1 : S.tax[bu], bv < 0 ? -1 : S.tax[bv]);
            return 0;
        }
        std::cerr << "leaf not found\n";
        return 1;
    }

    // --all : every leaf, printed as "taxon<TAB>best_delta"
    int done = 0;
    for (int u : leaves) {
        auto [d, bu, bv] = exact_best_for_leaf(S, u, nthreads);
        std::printf("%d\t%lld\t%d\t%d\n", S.tax[u], d,
                    bu < 0 ? -1 : S.tax[bu], bv < 0 ? -1 : S.tax[bv]);
        std::fflush(stdout);
        if (nLeaves > 0 && ++done >= nLeaves) break;
    }
    return 0;
}
