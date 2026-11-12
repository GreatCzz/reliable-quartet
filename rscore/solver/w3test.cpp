// spr: SPR (subtree prune-regraft) local search for the reliable-quartet objective.
//
// 与 solver.cpp 共用打分引擎（evalF / Wnode / W 缓存），区别在于邻域：
//   NNI: 每个内部节点试 2 种兄弟交换 -> 只影响 2 个节点的 W
//   SPR: 摘下一个子树 x，重挂到任意边 (u,v) -> 影响"两条根路径"上的节点 W
//
// 关键性质：SPR 只改变 x 到根路径（失去 x）和重挂点到根路径（获得 x）上节点的
// 子树/补集；两条路径的公共祖先"先失去再获得"净不变 -> 不用重算。
//
// 用法:
//   ./spr -s tree -i genes -o out.tree [--iters 2000] [--rounds 100] [--seed S] [-t 256]
//
// Build: g++ -std=c++20 -O3 -march=native -fopenmp -I ../aster2/src spr.cpp -o spr

#include "newick.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <omp.h>

namespace {

using std::size_t;
using std::string;
using std::vector;
using u64 = std::uint64_t;
using i64 = long long;

// ---------------------------------------------------------------- gene trees
struct GeneTree {
	vector<int> start, kids;
	vector<char> solid;
	vector<int> taxon;
	int root = 0, n = 0;
};

GeneTree buildGene(newick::Newick::Node const* const root) {
	GeneTree g;
	vector<newick::Newick::Node const*> order, st{ root };
	while (!st.empty()) {
		auto* u = st.back(); st.pop_back();
		order.push_back(u);
		for (auto it = u->children.rbegin(); it != u->children.rend(); ++it) st.push_back(*it);
	}
	std::unordered_map<newick::Newick::Node const*, int> id;
	id.reserve(order.size() * 2);
	for (size_t i = 0; i < order.size(); ++i) id[order[i]] = (int)i;
	g.n = (int)order.size();
	g.start.assign(g.n + 1, 0);
	g.solid.assign(g.n, 1);
	g.taxon.assign(g.n, -1);
	for (size_t i = 0; i < order.size(); ++i)
		g.start[i + 1] = g.start[i] + (int)order[i]->children.size();
	g.kids.reserve(g.start[g.n]);
	for (size_t i = 0; i < order.size(); ++i) {
		auto* u = order[i];
		if (u->children.empty()) g.taxon[i] = (int)common::taxonName2ID[u->name];
		else for (auto* c : u->children) {
			int const cid = id[c];
			g.kids.push_back(cid);
			g.solid[cid] = (c->name == "D") ? 0 : 1;
		}
	}
	g.root = 0;
	return g;
}

// --------------------------------------------------------------- F evaluation
struct Scratch {
	vector<u64> t1, t2;
	vector<i64> c1, a1, m1, G1, H1, c2, a2, m2, G2, H2;
	void resize(int n) {
		c1.resize(n); a1.resize(n); m1.resize(n); G1.resize(n); H1.resize(n);
		c2.resize(n); a2.resize(n); m2.resize(n); G2.resize(n); H2.resize(n);
	}
};

inline u64 evalF(GeneTree const& g, u64 const* const bl, u64 const* const br, Scratch& s) {
	int const n = g.n;
	i64* const c1 = s.c1.data(); i64* const a1 = s.a1.data(); i64* const m1 = s.m1.data();
	i64* const c2 = s.c2.data(); i64* const a2 = s.a2.data(); i64* const m2 = s.m2.data();
	for (int i = n - 1; i >= 0; --i) {
		int const t = g.taxon[i];
		if (t >= 0) {
			c1[i] = ((bl[(size_t)t >> 6] >> (t & 63)) & 1ULL) ? 1 : 0;
			c2[i] = ((br[(size_t)t >> 6] >> (t & 63)) & 1ULL) ? 1 : 0;
			a1[i] = 0; m1[i] = 0; a2[i] = 0; m2[i] = 0;
		}
		else {
			i64 C1 = 0, X1 = 0, M1 = 0, C2 = 0, X2 = 0, M2 = 0;
			int const e = g.start[i + 1];
			for (int k = g.start[i]; k < e; ++k) {
				int const ch = g.kids[k];
				M1 += m1[ch]; M2 += m2[ch];
				if (!g.solid[ch]) continue;
				i64 const v1 = c1[ch], v2 = c2[ch];
				X1 += C1 * v1; C1 += v1;
				X2 += C2 * v2; C2 += v2;
			}
			c1[i] = C1; a1[i] = X1; m1[i] = X1 + M1;
			c2[i] = C2; a2[i] = X2; m2[i] = X2 + M2;
		}
	}
	i64 const sp = m1[g.root], sq = m2[g.root];
	if (sp == 0 || sq == 0) return 0;
	i64* const G1 = s.G1.data(); i64* const H1 = s.H1.data();
	i64* const G2 = s.G2.data(); i64* const H2 = s.H2.data();
	G1[g.root] = 0; H1[g.root] = 0; G2[g.root] = 0; H2[g.root] = 0;
	u64 corr = (u64)a1[g.root] * (u64)a2[g.root];
	for (int i = 0; i < n; ++i) {
		int const e = g.start[i + 1];
		for (int k = g.start[i]; k < e; ++k) {
			int const ch = g.kids[k];
			if (g.solid[ch]) {
				G1[ch] = (c1[i] - c1[ch]) + G1[i];
				G2[ch] = (c2[i] - c2[ch]) + G2[i];
			}
			else { G1[ch] = 0; G2[ch] = 0; }
			H1[ch] = c1[ch] * G1[ch];
			H2[ch] = c2[ch] * G2[ch];
			corr += (u64)a1[ch] * (u64)a2[ch] + (u64)a2[ch] * (u64)H1[ch] + (u64)a1[ch] * (u64)H2[ch];
		}
	}
	u64 const total = (u64)sp * (u64)sq;
	return total > corr ? total - corr : 0;
}

inline i64 Wnode(GeneTree const& g, u64 const* A, u64 const* B, u64 const* C,
	int sA, int sB, int sC, size_t W, Scratch& sc) {
	i64 const fAB = (sA >= 2 && sB >= 2) ? (i64)evalF(g, A, B, sc) : 0;
	i64 const fAC = (sA >= 2 && sC >= 2) ? (i64)evalF(g, A, C, sc) : 0;
	i64 const fBC = (sB >= 2 && sC >= 2) ? (i64)evalF(g, B, C, sc) : 0;
	i64 gA = 0, gB = 0, gC = 0;
	u64* const t = sc.t2.data();
	if (sA >= 2 && sB + sC >= 2) {
		for (size_t w = 0; w < W; ++w) t[w] = ~A[w];
		gA = (i64)evalF(g, A, t, sc);
	}
	if (sB >= 2 && sA + sC >= 2) {
		for (size_t w = 0; w < W; ++w) t[w] = ~B[w];
		gB = (i64)evalF(g, B, t, sc);
	}
	if (sC >= 2 && sA + sB >= 2) {
		for (size_t w = 0; w < W; ++w) t[w] = ~C[w];
		gC = (i64)evalF(g, C, t, sc);
	}
	return gA + gB + gC - 2 * (fAB + fAC + fBC);
}

// --------------------------------------------------------------- species tree
struct SpeciesTree {
	int n = 0, root = 0;
	vector<int> par, c0, c1, taxon;
	vector<int> internal;
	vector<u64> bits;
	vector<int> sz;
	size_t W = 0;
};

SpeciesTree buildSpecies(newick::Newick::Node const* const root, size_t nTaxa) {
	SpeciesTree t;
	vector<newick::Newick::Node const*> order, st{ root };
	while (!st.empty()) {
		auto* u = st.back(); st.pop_back();
		order.push_back(u);
		for (auto it = u->children.rbegin(); it != u->children.rend(); ++it) st.push_back(*it);
	}
	std::unordered_map<newick::Newick::Node const*, int> id;
	id.reserve(order.size() * 2);
	for (size_t i = 0; i < order.size(); ++i) id[order[i]] = (int)i;
	t.n = (int)order.size();
	t.par.assign(t.n, -1); t.c0.assign(t.n, -1); t.c1.assign(t.n, -1); t.taxon.assign(t.n, -1);
	t.root = 0;
	t.W = (nTaxa + 63) / 64;
	t.bits.assign((size_t)t.n * t.W, 0);
	t.sz.assign(t.n, 0);
	for (size_t i = 0; i < order.size(); ++i) {
		auto* u = order[i];
		if (u->children.empty()) t.taxon[i] = (int)common::taxonName2ID[u->name];
		else {
			if (u->children.size() != 2) {
				std::cerr << "Error: species tree must be binary" << std::endl;
				std::exit(1);
			}
			t.c0[i] = id[u->children[0]];
			t.c1[i] = id[u->children[1]];
			t.par[t.c0[i]] = (int)i;
			t.par[t.c1[i]] = (int)i;
			t.internal.push_back((int)i);
		}
	}
	for (int i = t.n - 1; i >= 0; --i) {
		u64* const dst = t.bits.data() + (size_t)i * t.W;
		if (t.taxon[i] >= 0) {
			size_t const x = (size_t)t.taxon[i];
			dst[x >> 6] |= 1ULL << (x & 63);
			t.sz[i] = 1;
		}
		else {
			u64 const* const a = t.bits.data() + (size_t)t.c0[i] * t.W;
			u64 const* const b = t.bits.data() + (size_t)t.c1[i] * t.W;
			for (size_t w = 0; w < t.W; ++w) dst[w] = a[w] | b[w];
			t.sz[i] = t.sz[t.c0[i]] + t.sz[t.c1[i]];
		}
	}
	return t;
}

// 调试用：检查 par 链是否有环、根是否唯一
bool validateTree(SpeciesTree const& t) {
	for (int u = 0; u < t.n; ++u) {
		int steps = 0;
		for (int z = u; z >= 0; z = t.par[z]) {
			if (++steps > t.n) {
				std::cerr << "  [valid] cycle from " << u << std::endl;
				return false;
			}
		}
		if (t.par[u] >= 0 && t.c0[t.par[u]] != u && t.c1[t.par[u]] != u) {
			std::cerr << "  [valid] node " << u << " par=" << t.par[u]
				<< " but children are " << t.c0[t.par[u]] << "," << t.c1[t.par[u]] << std::endl;
			return false;
		}
	}
	return true;
}

string toNewick(SpeciesTree const& t, int u) {
	if (t.taxon[u] >= 0) return common::taxonName2ID[(size_t)t.taxon[u]];
	return "(" + toNewick(t, t.c0[u]) + "," + toNewick(t, t.c1[u]) + ")";
}


} // namespace

#include "w3.hpp"

int main(int argc, char** argv) {
	string spFile = argv[1], gFile = argv[2];
	int nG = (argc > 3) ? atoi(argv[3]) : 20;
	newick::Newick spTree;
	{
		std::ifstream fin(spFile);
		string line;
		while (std::getline(fin, line)) if (line.find_first_not_of(" \t\r\n") != string::npos) { spTree = newick::Newick(line); break; }
	}
	vector<GeneTree> genes;
	{
		std::ifstream fin(gFile);
		string line; int cnt = 0;
		while (std::getline(fin, line) && cnt < nG) {
			if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
			newick::Newick nw(line);
			genes.push_back(buildGene(nw.root()));
			++cnt;
		}
	}
	size_t const nTaxa = common::taxonName2ID.nTaxa();
	SpeciesTree S = buildSpecies(spTree.root(), nTaxa);
	int const T = (int)genes.size();
	vector<w3::GeneView> gv;
	for (auto& g : genes) {
		w3::GeneView v;
		v.start = g.start; v.kids = g.kids; v.solid = g.solid; v.taxon = g.taxon; v.root = g.root; v.n = g.n;
		gv.push_back(std::move(v));
	}
	int maxN = 1; for (auto& g : genes) maxN = std::max(maxN, g.n);
	Scratch sc; sc.resize(maxN); sc.t1.assign(S.W, 0); sc.t2.assign(S.W, 0);
	w3::Scratch3 s3; s3.resize(maxN); s3.t.assign(S.W, 0);
	long long bad = 0, checked = 0;
	double tOld = 0, tNew = 0;
	for (int k = 0; k < (int)S.internal.size(); ++k) {
		int const v = S.internal[(size_t)k];
		int const a = S.c0[v], b = S.c1[v];
		u64 const* const bA = S.bits.data() + (size_t)a * S.W;
		u64 const* const bB = S.bits.data() + (size_t)b * S.W;
		for (size_t w = 0; w < S.W; ++w) sc.t1[w] = ~(bA[w] | bB[w]);
		int const sA = S.sz[a], sB = S.sz[b], sC = (int)nTaxa - sA - sB;
		long long oldSum = 0, newSum = 0;
		auto t0 = std::chrono::steady_clock::now();
		long long fAB_ref_sum = 0, fAB_new_sum = 0, m_sum = 0;
		long long fABC_ref_sum = 0, fABC_new_sum = 0;
		i64 fABC_ref_sum_local = 0;
		for (int gi = 0; gi < T; ++gi) {
			oldSum += Wnode(genes[gi], bA, bB, sc.t1.data(), sA, sB, sC, S.W, sc);
			fAB_ref_sum += (sA>=2 && sB>=2) ? (i64)evalF(genes[gi], bA, bB, sc) : 0;
			{
				std::vector<u64> bR(S.W);
				for (size_t w = 0; w < S.W; ++w) bR[w] = bB[w] | sc.t1[w];
				int sR = (int)nTaxa - sA;
				i64 const v = (sA>=2 && sR>=2) ? (i64)evalF(genes[gi], bA, bR.data(), sc) : 0;
				fABC_ref_sum += v;
				fABC_ref_sum_local = v;
			}
			newSum += w3::Wnode3(gv[gi], bA, bB, sc.t1.data(), S.W, s3);
			fAB_new_sum += w3::g_dbg[0];
			fABC_new_sum += w3::g_dbg[3];
			if (gi == 0 && k == 3) {
				i64 fBC_ref = (sB>=2 && sC>=2) ? (i64)evalF(genes[gi], bB, sc.t1.data(), sc) : 0;
				i64 fAC_ref = (sA>=2 && sC>=2) ? (i64)evalF(genes[gi], bA, sc.t1.data(), sc) : 0;
				std::cerr << "    gene0: fBC_ref=" << fBC_ref << " fBC_new=" << w3::g_dbg[2]
					<< " | fAC_ref=" << fAC_ref << " fAC_new=" << w3::g_dbg[1]
					<< " | fA_BC_ref=" << fABC_ref_sum_local << "\n";
			}
		}
		auto t1 = std::chrono::steady_clock::now();
		// 单独再跑一遍 new（不计入 oldSum）
		long long dummy = 0;
		for (int gi = 0; gi < T; ++gi) dummy += w3::Wnode3(gv[gi], bA, bB, sc.t1.data(), S.W, s3);
		auto t2 = std::chrono::steady_clock::now();
		tOld += std::chrono::duration<double>(t1 - t0).count();
		tNew += std::chrono::duration<double>(t2 - t1).count();
		(void)dummy;
		++checked;
		if (k == 3 && oldSum != newSum) {
			std::cerr << "  node=" << v << " fAB_ref=" << fAB_ref_sum << " fAB_new=" << fAB_new_sum
				<< " | fA_BC_ref=" << fABC_ref_sum << " fA_BC_new=" << fABC_new_sum << "\n";
		}
		if (oldSum != newSum) {
			++bad;
			if (bad <= 5) std::cerr << "MISMATCH node=" << v << " old=" << oldSum << " new=" << newSum << "\n";
		}
	}
	std::cerr << "checked " << checked << " nodes x " << T << " genes, mismatches=" << bad << "\n";
	std::cerr << "Wnode(6x evalF) = " << tOld << " s | Wnode3(单遍) = " << tNew << " s | 加速 " << (tOld/tNew) << "x\n";
	return bad ? 1 : 0;
}
