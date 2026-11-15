// sector.hpp: 用"扇区聚合公式"一次算出"把叶子 x 插到每个位置"的 Δ（仅用前 G 棵树）
//
// 与 placedelta.hpp 等价，但直接工作在 spr 的 SpeciesTree/GeneTree 结构上。
#pragma once
#include <cstdint>
#include <unordered_map>
#include <functional>
#include <vector>
#include <algorithm>

namespace sector {

using i64 = long long;
using u64 = std::uint64_t;

// T: 物种树（需 par/c0/c1/taxon/bits/sz/root/W/n）；genes: 基因树 CSR；x: 目标叶子在 T 中的节点 id
// 返回 delta[u]（插到节点 u 上方，u != root）与 rootDelta（插到根之上）
template <typename PTree, typename GeneVec, typename ScratchT>
void computeDeltas(PTree const& T, GeneVec const& genes, int G, int x,
                   std::vector<i64>& delta, i64& rootDelta,
                   std::vector<int>& taxonOf, int nTaxa, ScratchT& sc) {
    int const n = T.n;
    delta.assign((size_t)n, 0);
    rootDelta = 0;
    if (x < 0 || T.taxon[x] < 0) return;
    int const taxonX = T.taxon[x];

    // 子树的叶子列表（用于 X 的子树和：直接按 DFS 序做）
    std::vector<int> ord;
    {
        std::vector<int> st{T.root};
        while (!st.empty()) { int u = st.back(); st.pop_back(); ord.push_back(u);
            if (T.c0[u] >= 0) { st.push_back(T.c0[u]); st.push_back(T.c1[u]); } }
    }

    std::vector<i64> sub((size_t)n, 0);
    i64 global = 0;

    for (int gi = 0; gi < G; ++gi) {
        auto const& g = genes[(size_t)gi];
        // ---- 1) 基因树分量（去掉 D 边）
        auto kidsOf = [&](int u, int& b, int& e) { b = g.start[(size_t)u]; e = g.start[(size_t)u + 1]; };
        std::vector<int> comp((size_t)g.n, -1);
        std::vector<int> compSize;
        std::function<void(int,int)> dfs = [&](int u, int cur) {
            int b, e; kidsOf(u, b, e);
            if (b == e) { comp[(size_t)u] = cur; return; }
            for (int k = b; k < e; ++k) {
                int const c = g.kids[(size_t)k];
                if (g.solid[(size_t)c]) dfs(c, cur);
                else { compSize.push_back(0); int nc = (int)compSize.size() - 1; dfs(c, nc); }
            }
        };
        compSize.push_back(0);
        dfs(g.root, 0);
        std::vector<int> cnt(compSize.size(), 0);
        for (int u = 0; u < g.n; ++u) { int b, e; kidsOf(u, b, e); if (b == e) cnt[(size_t)comp[(size_t)u]]++; }
        for (size_t k = 0; k < cnt.size(); ++k) compSize[k] = cnt[k];

        // ---- 2) 物种树叶 -> 分量；x 的分量
        std::vector<int> leafComp((size_t)n, -1);
        // 用 taxon -> 基因叶 的映射：直接遍历基因叶
        std::unordered_map<int,int> taxComp;
        taxComp.reserve(g.n * 2);
        for (int u = 0; u < g.n; ++u) { int b, e; kidsOf(u, b, e); if (b == e) taxComp[g.taxon[(size_t)u]] = comp[(size_t)u]; }
        for (int u = 0; u < n; ++u) if (T.taxon[u] >= 0) {
            auto it = taxComp.find(T.taxon[u]);
            if (it != taxComp.end()) leafComp[(size_t)u] = it->second;
        }
        auto itX = taxComp.find(taxonX);
        if (itX == taxComp.end()) continue;
        int const compX = itX->second;

        // ---- 3) X（子树和）与 degSum
        std::vector<i64> X((size_t)n, 0), degSum((size_t)n, 0);
        for (int u = 0; u < n; ++u) if (T.taxon[u] >= 0) {
            int const c = leafComp[(size_t)u];
            X[(size_t)u] = (c == compX) ? 1 : 0;
            degSum[(size_t)u] = (c >= 0) ? (compSize[(size_t)c] - 1) : 0;
            if (c >= 0 && c == compX && degSum[(size_t)u] > 0) --degSum[(size_t)u];
        }
        for (int i = (int)ord.size() - 1; i >= 0; --i) {
            int const u = ord[(size_t)i];
            if (T.taxon[u] >= 0) continue;
            X[(size_t)u] = X[(size_t)T.c0[u]] + X[(size_t)T.c1[u]];
            degSum[(size_t)u] = degSum[(size_t)T.c0[u]] + degSum[(size_t)T.c1[u]];
        }

        // ---- 4) pairsAt（LCA 为 m 的实线对数，小到大合并）
        std::vector<std::unordered_map<int,int>> maps((size_t)n);
        std::vector<i64> pairsAt((size_t)n, 0);
        for (int i = (int)ord.size() - 1; i >= 0; --i) {
            int const u = ord[(size_t)i];
            if (T.taxon[u] >= 0) {
                if (leafComp[(size_t)u] >= 0) maps[(size_t)u][leafComp[(size_t)u]] = 1;
                continue;
            }
            int a = T.c0[u], b = T.c1[u];
            if (maps[(size_t)a].size() < maps[(size_t)b].size()) std::swap(a, b);
            i64 cross = 0;
            for (auto const& kv : maps[(size_t)b]) {
                auto it = maps[(size_t)a].find(kv.first);
                if (it != maps[(size_t)a].end()) cross += (i64)it->second * kv.second;
            }
            for (auto const& kv : maps[(size_t)b]) maps[(size_t)a][kv.first] += kv.second;
            maps[(size_t)u] = std::move(maps[(size_t)a]);
            pairsAt[(size_t)u] = cross;
        }
        // ---- 5) pairsIn（子树内 pairsAt 之和）
        std::vector<i64> pairsIn((size_t)n, 0);
        for (int i = (int)ord.size() - 1; i >= 0; --i) {
            int const u = ord[(size_t)i];
            if (T.taxon[u] >= 0) { pairsIn[(size_t)u] = 0; continue; }
            pairsIn[(size_t)u] = pairsIn[(size_t)T.c0[u]] + pairsIn[(size_t)T.c1[u]] + pairsAt[(size_t)u];
        }
        // ---- 6) 每个内部节点的贡献（range add）
        for (int m = 0; m < n; ++m) {
            if (T.taxon[m] >= 0 || m == T.root) continue;
            int const c1 = T.c0[m], c2 = T.c1[m];
            i64 const X1 = X[(size_t)c1], X2 = X[(size_t)c2], X3 = X[(size_t)T.root] - X1 - X2;
            i64 const P12 = pairsAt[(size_t)m];
            i64 const P13 = degSum[(size_t)c1] - 2 * pairsIn[(size_t)c1] - P12;
            i64 const P23 = degSum[(size_t)c2] - 2 * pairsIn[(size_t)c2] - P12;
            sub[(size_t)c1] += X1 * P23;
            sub[(size_t)c2] += X2 * P13;
            global += X3 * P12;
            sub[(size_t)c1] -= X3 * P12;
            sub[(size_t)c2] -= X3 * P12;
        }
        (void)sc;
    }

    // ---- 7) 自上而下累加
    std::vector<i64> acc((size_t)n, 0);
    for (int i = 0; i < (int)ord.size(); ++i) {
        int const u = ord[(size_t)i];
        acc[(size_t)u] = (u == T.root) ? sub[(size_t)u]
                                       : acc[(size_t)T.par[u]] + sub[(size_t)u];
        delta[(size_t)u] = global + acc[(size_t)u];
    }
    rootDelta = global + acc[(size_t)T.root];
}

} // namespace sector
