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
    vector<string> name;
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


// 计算所有插入位置的 Δ（含"根的上一格"位置 rootPos）
// delta[u]：把 x 插到节点 u 上方（u != root）时的 Δ
// rootDelta：把 x 插到根之上（新建根）时的 Δ
long long computeAllDeltas(PTree const& T, vector<string> const& geneLines, int G,
                           string const& xname, vector<long long>& delta, long long& rootDelta) {
    int const nSp = (int)T.par.size();
    delta.assign(nSp, 0);
    rootDelta = 0;
    vector<string> nodeName(nSp);
    for (int u = 0; u < nSp; ++u) if (T.taxon[u] >= 0) nodeName[u] = T.name[(size_t)T.taxon[u]];
    // 子树叶集
    vector<vector<int>> subLeaves(nSp);
    vector<int> ord;
    { vector<int> st{T.root}; while (!st.empty()) { int u = st.back(); st.pop_back(); ord.push_back(u); if (T.c0[u] >= 0) { st.push_back(T.c0[u]); st.push_back(T.c1[u]); } } }
    for (int idx = (int)ord.size() - 1; idx >= 0; --idx) {
        int u = ord[idx];
        if (T.taxon[u] >= 0) subLeaves[u] = {u};
        else { subLeaves[u] = subLeaves[T.c0[u]]; subLeaves[u].insert(subLeaves[u].end(), subLeaves[T.c1[u]].begin(), subLeaves[T.c1[u]].end()); }
    }
    std::map<string,int> nameOfNode;
    for (int u = 0; u < nSp; ++u) if (T.taxon[u] >= 0) nameOfNode[nodeName[u]] = u;

    vector<long long> total(nSp, 0);
    for (int gi = 0; gi < G; ++gi) {
        GTree g;
        { int i = 0; while (i < (int)geneLines[gi].size() && geneLines[gi][i] != '(') ++i; if (i >= (int)geneLines[gi].size()) continue; g.root = parseNode(geneLines[gi], i, g); }
        int n = (int)g.ch.size();
        vector<char> isLeaf(n, 0);
        for (int u = 0; u < n; ++u) isLeaf[u] = g.ch[u].empty();
        vector<int> comp(n, -1), compSize;
        function<void(int,int)> dfs = [&](int u, int cur) {
            if (isLeaf[u]) { comp[u] = cur; return; }
            for (int c : g.ch[u]) {
                if (g.lab[c] == "D") { compSize.push_back(0); int nc = (int)compSize.size() - 1; dfs(c, nc); }
                else dfs(c, cur);
            }
        };
        compSize.push_back(0); dfs(g.root, 0);
        vector<int> cnt(compSize.size(), 0);
        for (int u = 0; u < n; ++u) if (isLeaf[u]) cnt[comp[u]]++;
        for (size_t k = 0; k < cnt.size(); ++k) compSize[k] = cnt[k];
        map<string,int> geneComp;
        for (int u = 0; u < n; ++u) if (isLeaf[u]) geneComp[g.lab[u]] = comp[u];
        auto itX = geneComp.find(xname);
        if (itX == geneComp.end()) continue;
        int const compX = itX->second;
        vector<int> leafComp(nSp, -1);
        for (int u = 0; u < nSp; ++u) if (T.taxon[u] >= 0) {
            auto it = geneComp.find(nodeName[u]);
            if (it != geneComp.end()) leafComp[u] = it->second;
        }
        vector<long long> X(nSp, 0), degSum(nSp, 0);
        for (int u = 0; u < nSp; ++u) if (T.taxon[u] >= 0) {
            X[u] = (leafComp[u] == compX) ? 1 : 0;
            degSum[u] = (leafComp[u] >= 0) ? (compSize[leafComp[u]] - 1) : 0;
            if (leafComp[u] >= 0 && leafComp[u] == compX && degSum[u] > 0) --degSum[u];
        }
        for (int idx = (int)ord.size() - 1; idx >= 0; --idx) {
            int u = ord[idx];
            if (T.taxon[u] >= 0) continue;
            X[u] = X[T.c0[u]] + X[T.c1[u]];
            degSum[u] = degSum[T.c0[u]] + degSum[T.c1[u]];
        }
        vector<map<int,int>> cntMap(nSp);
        vector<long long> pairsAt(nSp, 0);
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
        vector<long long> pairsIn(nSp, 0);
        for (int idx = (int)ord.size() - 1; idx >= 0; --idx) {
            int u = ord[idx];
            if (T.taxon[u] >= 0) { pairsIn[u] = 0; continue; }
            pairsIn[u] = pairsIn[T.c0[u]] + pairsIn[T.c1[u]] + pairsAt[u];
        }
        vector<long long> sub(nSp, 0);
        long long global = 0;
        for (int m = 0; m < nSp; ++m) {
            if (T.taxon[m] >= 0 || m == T.root) continue;
            int c1 = T.c0[m], c2 = T.c1[m];
            long long X1 = X[c1], X2 = X[c2], X3 = X[T.root] - X1 - X2;
            long long P12 = pairsAt[m];
            long long P13 = degSum[c1] - 2 * pairsIn[c1] - P12;
            long long P23 = degSum[c2] - 2 * pairsIn[c2] - P12;
            sub[c1] += X1 * P23;
            sub[c2] += X2 * P13;
            global += X3 * P12;
            sub[c1] -= X3 * P12;
            sub[c2] -= X3 * P12;
        }
        vector<long long> acc(nSp, 0);
        for (int idx = 0; idx < (int)ord.size(); ++idx) {
            int u = ord[idx];
            if (u == T.root) acc[u] = sub[u];
            else acc[u] = acc[T.par[u]] + sub[u];
            total[u] += global + acc[u];
        }
        rootDelta += global + acc[T.root];
    }
    delta = total;
    return rootDelta;
}
