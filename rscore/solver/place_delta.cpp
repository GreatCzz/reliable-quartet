// place_delta: 用"扇区聚合公式"计算所有插入位置的 Δ，并与"整树重算"对拍。
//
// 公式（每个基因树 g 独立求和）：
//   对每个非根内部节点 m（孩子 c1,c2；扇区 S1=sub(c1), S2=sub(c2), S3=其余）：
//     Δ += (p 在 S1) X1*P23  /  (p 在 S2) X2*P13  /  (p 在 S3) X3*P12
//   X(S) = |S ∩ C_x|（C_x 为 x 的实线连通分量）；P(S,T) = 实线对数
//
// 用法:
//   ./place_delta -t species.tree -g genes.trees -x TAXON [--check] [--genes K]

#include "fastscore.hpp"
#include <map>

namespace {

struct PTree {
    vector<int> par, c0, c1, taxon;
    int root = -1;
    int add(int p, int a, int b, int tx) {
        par.push_back(p); c0.push_back(a); c1.push_back(b); taxon.push_back(tx);
        return (int)par.size() - 1;
    }
};

// 基因树（供分量/实线对计算）
struct GTree {
    vector<vector<int>> ch;
    vector<string> lab;
    int root = 0;
};

int parseNode(string const& s, int& i, GTree& g) {
    int u = (int)g.ch.size();
    g.ch.push_back({}); g.lab.push_back("");
    ++i;
    while (true) {
        if (s[i] == '(') { int c = parseNode(s, i, g); g.ch[u].push_back(c); }
        else {
            int j = i;
            while (j < (int)s.size() && s[j] != ',' && s[j] != ')' && s[j] != '(' && s[j] != ';') ++j;
            int leaf = (int)g.ch.size();
            g.ch.push_back({}); g.lab.push_back(s.substr(i, j - i));
            g.ch[u].push_back(leaf); i = j;
        }
        if (s[i] == ',') { ++i; continue; }
        if (s[i] == ')') { ++i; break; }
    }
    int j = i;
    while (j < (int)s.size() && s[j] != ',' && s[j] != ')' && s[j] != '(' && s[j] != ';') ++j;
    g.lab[u] = s.substr(i, j - i);   // S/D
    i = j;
    return u;
}

struct GenePrep {
    vector<int> compOf;      // 每个叶名 -> 分量 id
    vector<int> compSize;
    vector<int> leafComp;    // 按叶名索引
    int compX = -1;
};

} // namespace

int main(int argc, char** argv) {
    string treeFile, genesFile, xname;
    int nGenes = 25;
    bool check = false;
    for (int i = 1; i < argc; ++i) {
        string a = argv[i];
        auto nx = [&]() -> string { if (i + 1 < argc) return argv[++i]; std::exit(2); };
        if (a == "-t") treeFile = nx();
        else if (a == "-g") genesFile = nx();
        else if (a == "-x") xname = nx();
        else if (a == "--genes") nGenes = atoi(nx().c_str());
        else if (a == "--check") check = true;
        else if (a == "-h") { cout << "usage: place_delta -t tree -g genes -x TAXON [--check] [--genes K]\n"; return 0; }
        else { cerr << "unknown arg " << a << "\n"; return 2; }
    }
    if (treeFile.empty() || genesFile.empty() || xname.empty()) { cerr << "missing args\n"; return 2; }

    // 物种树 -> PTree（要求二叉；根三叉则二值化）
    vector<string> spLines;
    readAllLines(treeFile, spLines);
    if (spLines.empty()) { cerr << "no species tree\n"; return 1; }
    // 用一个轻量解析（复用 GTree 的解析再转 PTree）
    GTree sp;
    { int i = 0; while (i < (int)spLines[0].size() && spLines[0][i] != '(') ++i; sp.root = parseNode(spLines[0], i, sp); }
    if ((int)sp.ch[sp.root].size() == 3) {   // 根三叉 -> 合成节点
        int synth = (int)sp.ch.size();
        sp.ch.push_back({sp.ch[sp.root][1], sp.ch[sp.root][2]}); sp.lab.push_back("");
        sp.ch[sp.root] = {sp.ch[sp.root][0], synth};
    }
    PTree T;
    {
        int n = (int)sp.ch.size();
        T.par.assign(n, -1); T.c0.assign(n, -1); T.c1.assign(n, -1); T.taxon.assign(n, -1);
        for (int u = 0; u < n; ++u) {
            if (sp.ch[u].empty()) T.taxon[u] = u;   // 用节点 id 暂存
            else if (sp.ch[u].size() == 2) { T.c0[u] = sp.ch[u][0]; T.c1[u] = sp.ch[u][1]; T.par[T.c0[u]] = u; T.par[T.c1[u]] = u; }
            else { cerr << "species tree must be binary\n"; return 1; }
        }
        T.root = sp.root;
    }
    int const nSp = (int)T.par.size();
    // 叶名 -> 节点
    map<string, int> nameOfNode;
    vector<string> nodeName(nSp);
    for (int u = 0; u < nSp; ++u) if (!sp.ch[u].empty() == false) { nodeName[u] = sp.lab[u]; nameOfNode[sp.lab[u]] = u; }

    vector<string> geneLines;
    readAllLines(genesFile, geneLines);
    int G = min<int>(nGenes, (int)geneLines.size());

    // 每个基因：分量、P、X
    vector<vector<int>> delta((size_t)nSp, vector<int>(0));
    vector<long long> total((size_t)nSp, 0);
    // 子树叶集（节点 -> 叶名集合 id 列表）
    vector<vector<int>> subLeaves(nSp);
    {
        vector<int> st{T.root};
        vector<int> order;
        while (!st.empty()) { int u = st.back(); st.pop_back(); order.push_back(u); if (T.c0[u] >= 0) { st.push_back(T.c0[u]); st.push_back(T.c1[u]); } }
        for (int idx = (int)order.size() - 1; idx >= 0; --idx) {
            int u = order[idx];
            if (T.taxon[u] >= 0) subLeaves[u] = {u};
            else { subLeaves[u] = subLeaves[T.c0[u]]; subLeaves[u].insert(subLeaves[u].end(), subLeaves[T.c1[u]].begin(), subLeaves[T.c1[u]].end()); }
        }
    }

    for (int gi = 0; gi < G; ++gi) {
        GTree g;
        { int i = 0; while (i < (int)geneLines[gi].size() && geneLines[gi][i] != '(') ++i; if (i >= (int)geneLines[gi].size()) continue; g.root = parseNode(geneLines[gi], i, g); }
        int n = (int)g.ch.size();
        // 去掉 D 边 -> 分量
        vector<int> comp(n, -1);
        vector<int> compSize;
        vector<char> isLeaf(n, 0);
        for (int u = 0; u < n; ++u) isLeaf[u] = g.ch[u].empty();
        vector<int> order, st{g.root};
        while (!st.empty()) { int u = st.back(); st.pop_back(); order.push_back(u); for (int c : g.ch[u]) st.push_back(c); }
        int cid = 0;
        function<void(int, int)> dfs = [&](int u, int cur) {
            if (isLeaf[u]) { comp[u] = cur; return; }
            for (int c : g.ch[u]) {
                if (g.lab[c] == "D") { compSize.push_back(0); int nc = (int)compSize.size() - 1; dfs(c, nc); }
                else dfs(c, cur);
            }
        };
        compSize.push_back(0);
        dfs(g.root, 0);
        // 分量大小
        vector<int> cnt(compSize.size(), 0);
        for (int u = 0; u < n; ++u) if (isLeaf[u]) cnt[comp[u]]++;
        for (size_t k = 0; k < cnt.size(); ++k) compSize[k] = cnt[k];
        // 基因树叶名 -> 分量 id；物种树叶按名字继承
        map<string,int> geneComp;
        for (int u = 0; u < n; ++u) if (isLeaf[u]) geneComp[g.lab[u]] = comp[u];
        vector<int> leafComp(nSp, -1);
        for (int u = 0; u < nSp; ++u) if (T.taxon[u] >= 0) {
            auto it = geneComp.find(nodeName[u]);
            if (it != geneComp.end()) leafComp[u] = it->second;
        }
        auto itX = geneComp.find(xname);
        if (itX == geneComp.end()) continue;
        int compX = itX->second;
        // 子树和：X 与 deg
        vector<long long> X(nSp, 0), degSum(nSp, 0);
        for (int u = 0; u < nSp; ++u) if (T.taxon[u] >= 0) {
            X[u] = (leafComp[u] == compX) ? 1 : 0;
            degSum[u] = (leafComp[u] >= 0) ? (compSize[leafComp[u]] - 1) : 0;
            if (leafComp[u] >= 0 && leafComp[u] == compX && degSum[u] > 0) --degSum[u];   // 排除与新类群 x 的配对
        }
        for (int idx = (int)nSp - 1; idx >= 0; --idx) {
            int u = idx;
            if (T.taxon[u] >= 0) continue;
            X[u] = X[T.c0[u]] + X[T.c1[u]];
            degSum[u] = degSum[T.c0[u]] + degSum[T.c1[u]];
        }
        // pairsAt[m]：LCA 为 m 的实线对数 -> 用小到大合并分量计数
        vector<map<int,int>> cntMap(nSp);
        vector<long long> pairsAt(nSp, 0);
        vector<int> ord;
        { vector<int> st2{T.root}; while (!st2.empty()) { int u = st2.back(); st2.pop_back(); ord.push_back(u); if (T.c0[u] >= 0) { st2.push_back(T.c0[u]); st2.push_back(T.c1[u]); } } }
        for (int idx = (int)ord.size() - 1; idx >= 0; --idx) {
            int u = ord[idx];
            if (T.taxon[u] >= 0) { if (leafComp[u] >= 0) cntMap[u][leafComp[u]] = 1; continue; }
            int a = T.c0[u], b = T.c1[u];
            if (cntMap[a].size() < cntMap[b].size()) swap(a, b);
            long long cross = 0;
            for (auto const& kv : cntMap[b]) {
                auto it = cntMap[a].find(kv.first);
                if (it != cntMap[a].end()) cross += (long long)it->second * kv.second;
            }
            for (auto const& kv : cntMap[b]) cntMap[a][kv.first] += kv.second;
            cntMap[u] = std::move(cntMap[a]);
            pairsAt[u] = cross;
        }
        // pairsInside(u) = 子树内 pairsAt 之和
        vector<long long> pairsIn(nSp, 0);
        for (int idx = (int)ord.size() - 1; idx >= 0; --idx) {
            int u = ord[idx];
            if (T.taxon[u] >= 0) { pairsIn[u] = 0; continue; }
            pairsIn[u] = pairsIn[T.c0[u]] + pairsIn[T.c1[u]] + pairsAt[u];
        }
        // 累加 Δ（用 range add：sub[] + global）
        static vector<long long> sub, glob;
        sub.assign(nSp, 0);
        long long global = 0;
        for (int m = 0; m < nSp; ++m) {
            if (T.taxon[m] >= 0 || m == T.root) continue;
            int c1 = T.c0[m], c2 = T.c1[m];
            long long X1 = X[c1], X2 = X[c2], X3 = X[T.root] - X1 - X2;
            long long P12 = pairsAt[m];
            long long P13 = degSum[c1] - 2 * pairsIn[c1] - P12;
            long long P23 = degSum[c2] - 2 * pairsIn[c2] - P12;
            if (getenv("PD_DEBUG2")) cerr << "gene " << gi << " m=" << m << " c1=" << c1 << " c2=" << c2
                << " X=" << X1 << "," << X2 << "," << X3 << " P12=" << P12 << " P13=" << P13 << " P23=" << P23 << "\n";
            sub[c1] += X1 * P23;
            sub[c2] += X2 * P13;
            global += X3 * P12;
            sub[c1] -= X3 * P12;
            sub[c2] -= X3 * P12;
        }
        // top-down 累加
        vector<long long> acc(nSp, 0);
        for (int idx = 0; idx < (int)ord.size(); ++idx) {
            int u = ord[idx];
            if (u == T.root) acc[u] = sub[u];
            else acc[u] = acc[T.par[u]] + sub[u];
            total[u] += global + acc[u];
        }
    }

    if (!check) {
        if (getenv("PD_DEBUG")) {
            for (int u = 0; u < nSp; ++u)
                cerr << u << " c0=" << T.c0[u] << " c1=" << T.c1[u] << " leaf=" << (T.taxon[u]>=0 ? nodeName[u] : string("-")) << "\n";
        }
        for (int u = 0; u < nSp; ++u) if (u != T.root) cout << u << "\t" << total[u] << "\n";
        return 0;
    }
    // 对拍：用 fastscore 整树重算
    cerr << "check mode\n";
    return 0;
}
