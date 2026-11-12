// Wnode3: 三分划 W(A,B,C) 的"单次遍历"版（替代 6 次 evalF）
//
// 数学：W = F(A,B∪C)+F(B,A∪C)+F(C,A∪B) − 2(F(A,B)+F(A,C)+F(B,C))
// 其中 F(L,R) = m_L[root]*m_R[root] − Σ_u [a_L[u]*a_R[u] + a_R[u]*H_L[u] + a_L[u]*H_R[u]]
// 三个颜色一次算清：c_X, a_XY, m_*, G_X, H_X
#pragma once

#include <cstdint>
#include <vector>
#include <algorithm>

namespace w3 {

using i64 = long long;
using u64 = std::uint64_t;

struct GeneView {
    std::vector<int> start, kids;   // CSR
    std::vector<char> solid;        // 边 (parent->node) 是否实线
    std::vector<int> taxon;         // 叶子的类群 id，内部为 -1
    int root = 0, n = 0;
};

struct Scratch3 {
    // 3 个颜色
    std::vector<i64> c[3];          // 子树内可与 u 实线相连的该色叶子数
    std::vector<i64> a[3][3];       // LCA 为 u 的该色对（同色/异色）
    std::vector<i64> m[3][3];       // 该色并集内的实线成对总数
    std::vector<i64> G[3], H[3];
    std::vector<u64> t;
    void resize(int n) {
        for (int i = 0; i < 3; ++i) { c[i].resize(n); G[i].resize(n); H[i].resize(n); }
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) { a[i][j].resize(n); m[i][j].resize(n); }
    }
};

// color of a taxon: 0=A, 1=B, 2=C, -1=不在任何集合
inline i64 Wnode3(GeneView const& g, u64 const* A, u64 const* B, u64 const* C, size_t W, Scratch3& s) {
    int const n = g.n;
    // ---- pass 1: bottom-up
    for (int i = n - 1; i >= 0; --i) {
        int const t = g.taxon[i];
        if (t >= 0) {
            int col = -1;
            if ((A[(size_t)t >> 6] >> (t & 63)) & 1ULL) col = 0;
            else if ((B[(size_t)t >> 6] >> (t & 63)) & 1ULL) col = 1;
            else if ((C[(size_t)t >> 6] >> (t & 63)) & 1ULL) col = 2;
            for (int k = 0; k < 3; ++k) { s.c[k][i] = (k == col) ? 1 : 0; s.G[k][i] = 0; s.H[k][i] = 0; }
            for (int x = 0; x < 3; ++x) for (int y = 0; y < 3; ++y) { s.a[x][y][i] = 0; s.m[x][y][i] = 0; }
            continue;
        }
        i64 C_[3] = {0,0,0};
        for (int x = 0; x < 3; ++x) for (int y = 0; y < 3; ++y) s.a[x][y][i] = 0;
        int const e = g.start[i + 1];
        for (int k = g.start[i]; k < e; ++k) {
            int const ch = g.kids[k];
            if (!g.solid[ch]) continue;
            for (int x = 0; x < 3; ++x) {
                for (int y = x; y < 3; ++y) {
                    if (x == y) s.a[x][x][i] += C_[x] * s.c[x][ch];
                    else s.a[x][y][i] += C_[x] * s.c[y][ch] + C_[y] * s.c[x][ch];
                }
            }
            for (int x = 0; x < 3; ++x) C_[x] += s.c[x][ch];
        }
        for (int x = 0; x < 3; ++x) s.c[x][i] = C_[x];
        // m：各色并集
        for (int x = 0; x < 3; ++x) {
            i64 tot = s.a[x][x][i];
            for (int k = g.start[i]; k < e; ++k) tot += s.m[x][x][g.kids[k]];   // 所有孩子（与实线无关）
            s.m[x][x][i] = tot;
        }
        // 两色并集 m_xy = a_xx + a_yy + a_xy + Σ m_xy[ch]
        static int const pairs[3][2] = {{0,1},{0,2},{1,2}};
        for (auto& pr : pairs) {
            int x = pr[0], y = pr[1];
            i64 tot = s.a[x][x][i] + s.a[y][y][i] + s.a[x][y][i];
            for (int k = g.start[i]; k < e; ++k) tot += s.m[x][y][g.kids[k]];   // 所有孩子
            s.m[x][y][i] = tot;
        }
    }
    // ---- pass 2: top-down G/H
    for (int k = 0; k < 3; ++k) { s.G[k][g.root] = 0; s.H[k][g.root] = 0; }
    for (int i = 0; i < n; ++i) {
        int const e = g.start[i + 1];
        for (int k = g.start[i]; k < e; ++k) {
            int const ch = g.kids[k];
            for (int x = 0; x < 3; ++x) {
                if (g.solid[ch]) s.G[x][ch] = (s.c[x][i] - s.c[x][ch]) + s.G[x][i];
                else s.G[x][ch] = 0;
                s.H[x][ch] = s.c[x][ch] * s.G[x][ch];
            }
        }
    }
    // ---- 组装 6 个 F
    // 每个节点的贡献：a_L*a_R + a_R*H_L + a_L*H_R
    // 用会合的并集：a_R = a_xx + a_yy + a_xy ; H_R = H_x+H_y + c_x*G_y + c_y*G_x
    i64 corr[6] = {0,0,0,0,0,0};   // 顺序：AB, AC, BC, A|BC, B|AC, C|AB
    for (int u = 0; u < n; ++u) {
        i64 const aAA = s.a[0][0][u], aBB = s.a[1][1][u], aCC = s.a[2][2][u];
        i64 const aAB = s.a[0][1][u], aAC = s.a[0][2][u], aBC = s.a[1][2][u];
        i64 const HA = s.H[0][u], HB = s.H[1][u], HC = s.H[2][u];
        i64 const cA = s.c[0][u], cB = s.c[1][u], cC = s.c[2][u];
        i64 const GA = s.G[0][u], GB = s.G[1][u], GC = s.G[2][u];
        auto term = [&](i64 aL, i64 HL, i64 aR, i64 HR) -> i64 { return aL*aR + aR*HL + aL*HR; };
        corr[0] += term(aAA, HA, aBB, HB);
        corr[1] += term(aAA, HA, aCC, HC);
        corr[2] += term(aBB, HB, aCC, HC);
        // A | B∪C
        i64 const aBC_ = aBB + aCC + aBC;
        i64 const HBC = HB + HC + cB*GC + cC*GB;
        corr[3] += term(aAA, HA, aBC_, HBC);
        // B | A∪C
        i64 const aAC_ = aAA + aCC + aAC;
        i64 const HAC = HA + HC + cA*GC + cC*GA;
        corr[4] += term(aBB, HB, aAC_, HAC);
        // C | A∪B
        i64 const aAB_ = aAA + aBB + aAB;
        i64 const HAB = HA + HB + cA*GB + cB*GA;
        corr[5] += term(aCC, HC, aAB_, HAB);
    }
    i64 const mA = s.m[0][0][g.root], mB = s.m[1][1][g.root], mC = s.m[2][2][g.root];
    i64 const mAB = s.m[0][1][g.root], mAC = s.m[0][2][g.root], mBC = s.m[1][2][g.root];
    auto F = [](i64 mL, i64 mR, i64 cr) -> i64 { i64 t = mL * mR - cr; return t > 0 ? t : 0; };
    i64 const fAB = F(mA, mB, corr[0]);
    i64 const fAC = F(mA, mC, corr[1]);
    i64 const fBC = F(mB, mC, corr[2]);
    i64 const fA_BC = F(mA, mBC, corr[3]);
    i64 const fB_AC = F(mB, mAC, corr[4]);
    i64 const fC_AB = F(mC, mAB, corr[5]);
    return fA_BC + fB_AC + fC_AB - 2 * (fAB + fAC + fBC);
}

} // namespace w3
