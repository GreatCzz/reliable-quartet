// placer: 慢速但精确的"整树评估 + 打分放置"物种树构造器。
//
// 思路：逐个插入类群；每个候选位置都**实际插入并整树重新计分**（用导师快速 DP），取最优。
// 不做任何增量 Δ，避免深树/缺类群下的近似与 bug。
//
// 用法:
//   ./placer -i genes.trees -o out.tree [-t THREADS] [--genes K] [--cand C]
//            [--backbone B] [--order freq|rand] [--seed S]
//
//   --genes K    决策时使用的基因树数量（默认 20；最终输出会用全部基因树打分）
//   --cand C     每个类群的候选位置数（默认 6）
//   --backbone B 前 B 个类群用"全部位置"放置（默认 200）

#include "fastscore.hpp"

namespace {

struct PTree {
    vector<int> par, c0, c1, taxon;   // taxon >= 0 为叶
    int root = -1;
    int add(int p, int a, int b, int tx) {
        par.push_back(p); c0.push_back(a); c1.push_back(b); taxon.push_back(tx);
        return (int)par.size() - 1;
    }
};

// 基因树小解析器（仅用于算"邻居"）
struct GTree {
    vector<vector<int>> ch;
    vector<string> lab;
    int root = 0;
};

int parseGeneNode(string const& s, int& i, GTree& g) {
    int u = (int)g.ch.size();
    g.ch.push_back({}); g.lab.push_back("");
    ++i;  // '('
    while (true) {
        if (i < (int)s.size() && s[i] == '(') {
            int c = parseGeneNode(s, i, g);
            g.ch[u].push_back(c);
        } else {
            int j = i;
            while (j < (int)s.size() && s[j] != ',' && s[j] != ')' && s[j] != '(' && s[j] != ';') ++j;
            int leaf = (int)g.ch.size();
            g.ch.push_back({}); g.lab.push_back(s.substr(i, j - i));
            g.ch[u].push_back(leaf);
            i = j;
        }
        if (i < (int)s.size() && s[i] == ',') { ++i; continue; }
        if (i < (int)s.size() && s[i] == ')') { ++i; break; }
    }
    int j = i;
    while (j < (int)s.size() && s[j] != ',' && s[j] != ')' && s[j] != '(' && s[j] != ';') ++j;
    g.lab[u] = s.substr(i, j - i);   // 边标号 S/D
    i = j;
    return u;
}

// 计算每个 leaf 在给定基因树里的"最近邻居"（子树里另一个叶子的名字）
void collectPartners(string const& line, unordered_map<string, vector<string>>& partners, int capPerTaxon) {
    GTree g;
    int i = 0;
    while (i < (int)line.size() && line[i] != '(') ++i;
    if (i >= (int)line.size()) return;
    g.root = parseGeneNode(line, i, g);
    int n = (int)g.ch.size();
    vector<int> sz(n, 0);
    vector<int> firstLeaf(n, -1);
    function<void(int)> dfs = [&](int u) {
        if (g.ch[u].empty()) { sz[u] = 1; firstLeaf[u] = u; return; }
        for (int c : g.ch[u]) { dfs(c); sz[u] += sz[c]; if (firstLeaf[u] < 0) firstLeaf[u] = firstLeaf[c]; }
    };
    dfs(g.root);
    // 父指针
    vector<int> par(n, -1);
    for (int u = 0; u < n; ++u) for (int c : g.ch[u]) par[c] = u;
    for (int u = 0; u < n; ++u) {
        if (!g.ch[u].empty()) continue;
        string const& name = g.lab[u];
        vector<string>& lst = partners[name];
        if ((int)lst.size() >= capPerTaxon) continue;
        int cur = par[u];
        while (cur >= 0) {
            if (sz[cur] >= 2) {
                // 找子树里另一个叶子
                int other = -1;
                function<bool(int)> find = [&](int x) -> bool {
                    if (g.ch[x].empty()) { if (x != u) { other = x; return true; } return false; }
                    for (int c : g.ch[x]) if (find(c)) return true;
                    return false;
                };
                find(cur);
                if (other >= 0) {
                    string const& nm = g.lab[other];
                    bool dup = false;
                    for (auto const& z : lst) if (z == nm) { dup = true; break; }
                    if (!dup) lst.push_back(nm);
                }
                break;
            }
            cur = par[cur];
        }
    }
}

string toNewick(PTree const& T, int u, vector<string> const& names) {
    if (T.taxon[u] >= 0) return names[(size_t)T.taxon[u]];
    return "(" + toNewick(T, T.c0[u], names) + "," + toNewick(T, T.c1[u], names) + ")";
}

} // namespace

int main(int argc, char** argv) {
    string genesFile, outFile;
    int threads = 0, nDecide = 20, nCand = 6, backbone = 200;
    long seed = -1;
    bool randomOrder = false;
    for (int i = 1; i < argc; ++i) {
        string a = argv[i];
        auto nx = [&]() -> string { if (i + 1 < argc) return argv[++i]; std::exit(2); };
        if (a == "-i") genesFile = nx();
        else if (a == "-o") outFile = nx();
        else if (a == "-t") threads = atoi(nx().c_str());
        else if (a == "--genes") nDecide = atoi(nx().c_str());
        else if (a == "--cand") nCand = atoi(nx().c_str());
        else if (a == "--backbone") backbone = atoi(nx().c_str());
        else if (a == "--seed") seed = atol(nx().c_str());
        else if (a == "--order") { string v = nx(); randomOrder = (v == "rand"); }
        else if (a == "-h" || a == "--help") {
            cout << "usage: placer -i genes.trees -o out.tree [-t N] [--genes K] [--cand C] [--backbone B] [--order freq|rand] [--seed S]\n";
            return 0;
        }
        else { cerr << "Error: unknown argument " << a << "\n"; return 2; }
    }
    if (genesFile.empty()) { cerr << "Error: -i <gene trees> is required\n"; return 2; }

    vector<string> geneLines;
    if (!readAllLines(genesFile, geneLines)) { cerr << "Error: cannot open " << genesFile << "\n"; return 1; }
    if (geneLines.empty()) { cerr << "Error: no gene trees\n"; return 1; }

    // 类群表与频次
    unordered_map<string, int> name2id;
    vector<string> names;
    vector<long long> freq;
    for (auto const& line : geneLines) {
        GTree g;
        int i = 0; while (i < (int)line.size() && line[i] != '(') ++i;
        if (i >= (int)line.size()) continue;
        g.root = parseGeneNode(line, i, g);
        for (int u = 0; u < (int)g.ch.size(); ++u) if (g.ch[u].empty()) {
            auto it = name2id.find(g.lab[u]);
            int id;
            if (it == name2id.end()) { id = (int)names.size(); name2id.emplace(g.lab[u], id); names.push_back(g.lab[u]); freq.push_back(0); }
            else id = it->second;
            freq[(size_t)id]++;
        }
    }
    int const nTaxa = (int)names.size();
    if (nTaxa < 4) { cerr << "Error: fewer than 4 taxa\n"; return 1; }
    cerr << "[placer] taxa=" << nTaxa << " genes=" << geneLines.size()
         << " decide=" << min<int>(nDecide, (int)geneLines.size()) << "\n";

    // 邻居表
    unordered_map<string, vector<string>> partnerNames;
    for (int gi = 0; gi < min<int>(50, (int)geneLines.size()); ++gi)
        collectPartners(geneLines[gi], partnerNames, 3);
    vector<vector<int>> partners((size_t)nTaxa);
    for (int t = 0; t < nTaxa; ++t) {
        auto it = partnerNames.find(names[(size_t)t]);
        if (it == partnerNames.end()) continue;
        for (auto const& nm : it->second) {
            auto jt = name2id.find(nm);
            if (jt != name2id.end() && jt->second != t) partners[(size_t)t].push_back(jt->second);
        }
    }

    vector<int> order((size_t)nTaxa);
    for (int i = 0; i < nTaxa; ++i) order[(size_t)i] = i;
    if (!randomOrder)
        sort(order.begin(), order.end(), [&](int a, int b) { return freq[(size_t)a] > freq[(size_t)b]; });
    if (seed >= 0) { std::mt19937_64 rng((uint64_t)seed); std::shuffle(order.begin(), order.end(), rng); }

    PTree T;
    spPool.clear();
    spRootBase = -1; spN = 0;

    auto applyInsert = [&](int p, int x) -> array<int,4> {   // {par0, slot, oldRoot, n0}
        int par0 = T.par[p];
        int oldRoot = T.root;
        int n0 = (int)T.par.size();
        int leaf = n0, w = n0 + 1;
        T.par.push_back(par0); T.c0.push_back(-1); T.c1.push_back(-1); T.taxon.push_back(x);
        T.par.push_back(par0); T.c0.push_back(p); T.c1.push_back(leaf); T.taxon.push_back(-1);
        int slot = -1;
        if (par0 >= 0) {
            if (T.c0[par0] == p) { T.c0[par0] = w; slot = 0; } else { T.c1[par0] = w; slot = 1; }
            T.par[p] = w;
        } else { T.par[p] = w; T.root = w; }

        spPool.push_back(Node()); spPool.back().label = names[(size_t)x];
        spPool.push_back(Node()); spPool.back().ch.clear(); spPool.back().ch.push_back(p); spPool.back().ch.push_back(leaf);
        if (par0 >= 0) {
            spPool[par0].ch.clear(); spPool[par0].ch.push_back(T.c0[par0]); spPool[par0].ch.push_back(T.c1[par0]);
        }
        spRootBase = T.root; spN = (int)spPool.size(); g_spVersion++;
        return {par0, slot, oldRoot, n0};
    };
    auto undoInsert = [&](array<int,4> const& u, int p) {
        int par0 = u[0], slot = u[1], oldRoot = u[2], n0 = u[3];
        T.par.resize(n0); T.c0.resize(n0); T.c1.resize(n0); T.taxon.resize(n0);
        T.par[p] = par0;
        if (par0 >= 0) { if (slot == 0) T.c0[par0] = p; else T.c1[par0] = p; }
        T.root = oldRoot;
        spPool.resize(n0);
        if (par0 >= 0) {
            spPool[par0].ch.clear(); spPool[par0].ch.push_back(T.c0[par0]); spPool[par0].ch.push_back(T.c1[par0]);
        }
        spRootBase = T.root; spN = (int)spPool.size(); g_spVersion++;
    };

    int const nThreads = threads > 0 ? threads : (int)std::thread::hardware_concurrency();
    long long evals = 0;
    auto t0 = chrono::steady_clock::now();

    for (int idx = 0; idx < nTaxa; ++idx) {
        int const x = order[(size_t)idx];
        if (T.par.empty()) {
            T.add(-1, -1, -1, x);
            spPool.push_back(Node()); spPool.back().label = names[(size_t)x];
            T.root = 0; spRootBase = 0; spN = 1; g_spVersion++;
            continue;
        }

        vector<int> cand;
        if (idx < backbone) {
            cand.resize(T.par.size());
            for (int v = 0; v < (int)T.par.size(); ++v) cand[(size_t)v] = v;
        } else {
            cand.push_back(T.root);
            for (int y : partners[(size_t)x]) {
                int ny = -1;
                for (int v = 0; v < (int)T.par.size(); ++v) if (T.taxon[v] == y) { ny = v; break; }
                if (ny < 0) continue;
                cand.push_back(ny);
                int p = T.par[ny];
                for (int k = 0; k < 3 && p >= 0; ++k) { cand.push_back(p); p = T.par[p]; }
            }
            sort(cand.begin(), cand.end());
            cand.erase(unique(cand.begin(), cand.end()), cand.end());
            if ((int)cand.size() > nCand) cand.resize(nCand);
        }

        long long bestScore = LLONG_MIN;
        int bestP = cand.empty() ? T.root : cand[0];
        for (int p : cand) {
            auto u = applyInsert(p, x);
            long long sc = scoreSpeciesWithPool(geneLines, nDecide, nThreads);
            ++evals;
            undoInsert(u, p);
            if (sc > bestScore) { bestScore = sc; bestP = p; }
        }
        applyInsert(bestP, x);   // 正式插入（不撤销）

        if ((idx + 1) % 200 == 0 || idx + 1 == nTaxa)
            cerr << "[placer] " << (idx + 1) << "/" << nTaxa << "  evals=" << evals
                 << "  (" << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - t0).count() << " s)\n";
    }

    string nw = toNewick(T, T.root, names) + ";";
    if (outFile.empty()) cout << nw << "\n";
    else { ofstream f(outFile); f << nw << "\n"; }
    cerr << "[placer] done: " << nTaxa << " taxa, evals=" << evals << "\n";
    return 0;
}
