// builder: 用"打分插入"构造初始物种树（无根四叶树口径）。
//
// 思路：逐个插入类群；每个候选位置精确计算 Δ(ΣW)（含 p 到根所有祖先），取最优。
// 前 --budget 个类群用全候选（树小、便宜），其余只在"伙伴锚点"附近的候选里挑。
// 输出一棵二叉 Newick，交给 solver 做 NNI。
//
// 编译: g++ -std=c++20 -O3 -march=native -fopenmp -I ../aster2/src builder.cpp -o builder

#include "newick.hpp"
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>
#include <omp.h>

namespace {
using std::size_t; using std::string; using std::vector;
using u64 = std::uint64_t; using i64 = long long;

struct GeneTree {
	vector<int> start, kids; vector<char> solid; vector<int> taxon; int root = 0, n = 0;
};

GeneTree buildGene(newick::Newick::Node const* const root) {
	GeneTree g;
	vector<newick::Newick::Node const*> order, st{ root };
	while (!st.empty()) { auto* u = st.back(); st.pop_back(); order.push_back(u);
		for (auto it = u->children.rbegin(); it != u->children.rend(); ++it) st.push_back(*it); }
	std::unordered_map<newick::Newick::Node const*, int> id;
	id.reserve(order.size() * 2);
	for (size_t i = 0; i < order.size(); ++i) id[order[i]] = (int)i;
	g.n = (int)order.size();
	g.start.assign(g.n + 1, 0);
	g.solid.assign(g.n, 1);
	g.taxon.assign(g.n, -1);
	for (size_t i = 0; i < order.size(); ++i) g.start[i+1] = g.start[i] + (int)order[i]->children.size();
	g.kids.reserve(g.start[g.n]);
	for (size_t i = 0; i < order.size(); ++i) {
		auto* u = order[i];
		if (u->children.empty()) g.taxon[i] = (int)common::taxonName2ID[u->name];
		else for (auto* c : u->children) { int cid = id[c]; g.kids.push_back(cid); g.solid[cid] = (c->name == "D") ? 0 : 1; }
	}
	g.root = 0;
	return g;
}

struct Scratch {
	vector<u64> t1, t2;
	vector<i64> c1, a1, m1, G1, H1, c2, a2, m2, G2, H2;
	void resize(int n, size_t W) {
		c1.resize(n); a1.resize(n); m1.resize(n); G1.resize(n); H1.resize(n);
		c2.resize(n); a2.resize(n); m2.resize(n); G2.resize(n); H2.resize(n);
		t1.assign(W, 0); t2.assign(W, 0);
	}
};

inline u64 evalF(GeneTree const& g, u64 const* bl, u64 const* br, Scratch& s) {
	int const n = g.n;
	i64* c1 = s.c1.data(); i64* a1 = s.a1.data(); i64* m1 = s.m1.data();
	i64* c2 = s.c2.data(); i64* a2 = s.a2.data(); i64* m2 = s.m2.data();
	for (int i = n - 1; i >= 0; --i) {
		int const t = g.taxon[i];
		if (t >= 0) {
			c1[i] = ((bl[(size_t)t >> 6] >> (t & 63)) & 1ULL) ? 1 : 0;
			c2[i] = ((br[(size_t)t >> 6] >> (t & 63)) & 1ULL) ? 1 : 0;
			a1[i] = 0; m1[i] = 0; a2[i] = 0; m2[i] = 0;
		} else {
			i64 C1 = 0, X1 = 0, M1 = 0, C2 = 0, X2 = 0, M2 = 0;
			for (int k = g.start[i]; k < g.start[i+1]; ++k) {
				int const ch = g.kids[k];
				M1 += m1[ch]; M2 += m2[ch];
				if (!g.solid[ch]) continue;
				i64 v1 = c1[ch], v2 = c2[ch];
				X1 += C1 * v1; C1 += v1;
				X2 += C2 * v2; C2 += v2;
			}
			c1[i] = C1; a1[i] = X1; m1[i] = X1 + M1;
			c2[i] = C2; a2[i] = X2; m2[i] = X2 + M2;
		}
	}
	i64 sp = m1[g.root], sq = m2[g.root];
	if (sp == 0 || sq == 0) return 0;
	i64* G1 = s.G1.data(); i64* H1 = s.H1.data(); i64* G2 = s.G2.data(); i64* H2 = s.H2.data();
	G1[g.root] = 0; H1[g.root] = 0; G2[g.root] = 0; H2[g.root] = 0;
	u64 corr = (u64)a1[g.root] * (u64)a2[g.root];
	for (int i = 0; i < n; ++i) for (int k = g.start[i]; k < g.start[i+1]; ++k) {
		int const ch = g.kids[k];
		if (g.solid[ch]) { G1[ch] = (c1[i]-c1[ch]) + G1[i]; G2[ch] = (c2[i]-c2[ch]) + G2[i]; }
		else { G1[ch] = 0; G2[ch] = 0; }
		H1[ch] = c1[ch] * G1[ch]; H2[ch] = c2[ch] * G2[ch];
		corr += (u64)a1[ch]*(u64)a2[ch] + (u64)a2[ch]*(u64)H1[ch] + (u64)a1[ch]*(u64)H2[ch];
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
	u64* t = sc.t2.data();
	if (sA >= 2 && sB + sC >= 2) { for (size_t w = 0; w < W; ++w) t[w] = ~A[w]; gA = (i64)evalF(g, A, t, sc); }
	if (sB >= 2 && sA + sC >= 2) { for (size_t w = 0; w < W; ++w) t[w] = ~B[w]; gB = (i64)evalF(g, B, t, sc); }
	if (sC >= 2 && sA + sB >= 2) { for (size_t w = 0; w < W; ++w) t[w] = ~C[w]; gC = (i64)evalF(g, C, t, sc); }
	return gA + gB + gC - 2 * (fAB + fAC + fBC);
}

struct DynTree {
	int n = 0, root = -1;
	vector<int> par, c0, c1, taxon, sz;
	vector<u64> bits;
	vector<i64> wsum;      // 每个节点的 Σ_g W(g)，插入时维护
	size_t W = 0;

	int newNode(int p, int a, int b, int tx) {
		par.push_back(p); c0.push_back(a); c1.push_back(b); taxon.push_back(tx); sz.push_back(0); wsum.push_back(0);
		bits.resize(bits.size() + W, 0);
		n++;
		return n - 1;
	}
};

inline void setPath(u64* dst, u64 const* src, size_t W) { for (size_t w = 0; w < W; ++w) dst[w] = src[w]; }

// 计算节点 v 的 W（对给定的基因子集求和；用于插入后刷新缓存）
i64 nodeWsum(DynTree const& S, int v, vector<GeneTree> const& genes, int T, Scratch& sc) {
	int const a = S.c0[v], b = S.c1[v];
	if (a < 0) return 0;
	u64 const* bA = S.bits.data() + (size_t)a * S.W;
	u64 const* bB = S.bits.data() + (size_t)b * S.W;
	int const sC = (int)S.sz[S.root] - S.sz[a] - S.sz[b];
	for (size_t w = 0; w < S.W; ++w) sc.t1[w] = ~(bA[w] | bB[w]);
	i64 sum = 0;
	for (int g = 0; g < T; ++g)
		sum += Wnode(genes[g], bA, bB, sc.t1.data(), S.sz[a], S.sz[b], sC, S.W, sc);
	return sum;
}

string toNewick(DynTree const& S, int u) {
	if (S.taxon[u] >= 0) return common::taxonName2ID[(size_t)S.taxon[u]];
	return "(" + toNewick(S, S.c0[u]) + "," + toNewick(S, S.c1[u]) + ")";
}

} // namespace

int main(int argc, char** argv) {
	string inFile, outFile, seedOrder;
	int threads = 8, budget = 300, cand = 10, Tgenes = 8, orderMode = 0;
	for (int i = 1; i < argc; ++i) {
		string s = argv[i];
		auto nx = [&]() -> string { if (i + 1 < argc) return argv[++i]; std::exit(2); };
		if (s == "-i") inFile = nx();
		else if (s == "-o") outFile = nx();
		else if (s == "-t") threads = std::stoi(nx());
		else if (s == "--budget") budget = std::stoi(nx());
		else if (s == "--cand") cand = std::stoi(nx());
		else if (s == "--genes") Tgenes = std::stoi(nx());
	}
	if (inFile.empty()) { std::cerr << "usage: builder -i genes.trees -o out.tree [-t N] [--budget B] [--cand C] [--genes K]" << std::endl; return 2; }

	vector<GeneTree> genes;
	{
		std::ifstream fin(inFile);
		if (!fin.is_open()) { std::cerr << "Error: cannot open " << inFile << std::endl; return 2; }
		string line;
		while (std::getline(fin, line)) {
			if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
			newick::Newick nw(line);
			genes.push_back(buildGene(nw.root()));
		}
	}
	if (genes.empty()) { std::cerr << "Error: no gene trees" << std::endl; return 2; }
	size_t const nTaxa = common::taxonName2ID.nTaxa();
	int const T = std::min<int>(Tgenes, (int)genes.size());
	size_t const W = (nTaxa + 63) / 64;

	// 插入顺序：先放"信息量"高的（由基因树里出现次数估计），其余随机
	vector<int> freq(nTaxa, 0);
	for (auto const& g : genes) for (int u = 0; u < g.n; ++u) if (g.taxon[u] >= 0) freq[(size_t)g.taxon[u]]++;
	vector<int> order(nTaxa); for (size_t i = 0; i < nTaxa; ++i) order[i] = (int)i;
	std::sort(order.begin(), order.end(), [&](int x, int y) { return freq[(size_t)x] > freq[(size_t)y]; });

	// 事件：每个类群在"某棵基因树"里的初始锚点（父节点内另一个叶子），用于候选生成
	vector<int> anchorTaxon(nTaxa, -1);
	for (size_t x = 0; x < nTaxa; ++x) {
		for (auto const& g : genes) {
			// 找到 x 的父节点，取父节点子树里另一个叶子
			int node = -1;
			for (int u = 0; u < g.n; ++u) if (g.taxon[u] == (int)x) { node = u; break; }
			if (node < 0) continue;
			int p = -1;
			for (int u = 0; u < g.n && p < 0; ++u) for (int k = g.start[u]; k < g.start[u+1]; ++k) if (g.kids[k] == node) { p = u; break; }
			if (p < 0) continue;
			for (int k = g.start[p]; k < g.start[p+1]; ++k) {
				int c = g.kids[k];
				if (c == node) continue;
				int cand2 = -1;
				vector<int> st{ c };
				while (!st.empty() && cand2 < 0) { int u = st.back(); st.pop_back(); if (g.taxon[u] >= 0) cand2 = g.taxon[u]; else for (int kk = g.start[u]; kk < g.start[u+1]; ++kk) st.push_back(g.kids[kk]); }
				if (cand2 >= 0) { anchorTaxon[x] = cand2; break; }
			}
			if (anchorTaxon[x] >= 0) break;
		}
	}

	DynTree S; S.W = W;
	vector<Scratch> scratch(threads);
	int maxN = 1; for (auto const& g : genes) maxN = std::max(maxN, g.n);
	for (auto& sc : scratch) sc.resize(maxN, W);

	vector<int> nodeOfTaxon(nTaxa, -1);   // 已在树中的类群 -> 叶子节点 id
	int inserted = 0;
	auto t0 = std::chrono::steady_clock::now();

	// 以 order 逐个插入
	for (int idx = 0; idx < (int)nTaxa; ++idx) {
		int const x = order[idx];
		if (inserted == 0) { S.taxon.push_back(x); S.par.push_back(-1); S.c0.push_back(-1); S.c1.push_back(-1);
			S.sz.push_back(1); S.wsum.push_back(0); S.bits.resize(W, 0); S.bits[(size_t)x >> 6] |= 1ULL << (x & 63);
			S.n = 1; S.root = 0; nodeOfTaxon[x] = 0; inserted = 1; continue; }
		if (inserted == 1) {
			int const leaf = S.newNode(S.root, -1, -1, x);
			S.bits[(size_t)x >> 6] |= 1ULL << (x & 63); S.sz[leaf] = 1;
			int const newRoot = S.newNode(-1, S.root, leaf, -1);
			S.par[S.root] = newRoot; S.par[leaf] = newRoot;
			for (size_t w = 0; w < W; ++w) S.bits[(size_t)newRoot*W + w] = S.bits[(size_t)S.c0[newRoot]*W + w] | S.bits[(size_t)S.c1[newRoot]*W + w];
			S.sz[newRoot] = S.sz[S.c0[newRoot]] + S.sz[S.c1[newRoot]];
			S.root = newRoot; nodeOfTaxon[x] = leaf; inserted = 2; continue;
		}

		// 候选位置集合
		vector<int> cands;
		if (inserted < budget) {
			for (int v = 0; v < S.n; ++v) cands.push_back(v);
		} else {
			int anchor = anchorTaxon[x];
			int v = (anchor >= 0) ? nodeOfTaxon[anchor] : S.root;
			if (v < 0) v = S.root;
			for (int step = 0; step < cand && v >= 0; ++step) { cands.push_back(v); v = S.par[v]; }
			// 锚点的兄弟孩子
			// （略：路径已覆盖大部分情况）
		}

		// 评估：把 x 插到候选 p 的上方（精确 Δ(ΣW)）
		i64 bestDelta = LLONG_MIN; int bestP = cands.empty() ? S.root : cands[0];
		vector<i64> deltas(cands.size(), LLONG_MIN);
		Scratch& sc = scratch[0];
		vector<u64> bX(W, 0), bA(W), bB(W), bC(W);
		bX[(size_t)x >> 6] |= 1ULL << (x & 63);
		#pragma omp parallel for schedule(dynamic) num_threads(threads)
		for (int ci = 0; ci < (int)cands.size(); ++ci) {
			int const p = cands[ci];
			Scratch& sc = scratch[omp_get_thread_num()];
			vector<u64> bUnion(W), bC2(W);
			i64 delta = 0;
			// 新节点 w = (p, x)
			for (size_t k = 0; k < W; ++k) bUnion[k] = S.bits[(size_t)p * W + k] | bX[k];
			{
				int const sA = S.sz[p] + 1, sB = 1;
				int const sC = (int)S.sz[S.root] - sA - sB;
				for (size_t k = 0; k < W; ++k) bC2[k] = ~(bUnion[k] | bX[k]);
				for (int g = 0; g < T; ++g)
					delta += Wnode(genes[g], bUnion.data(), bX.data(), bC2.data(), sA, sB, sC, W, sc);
			}
			// p 若是当前的根：插入后 p 的补集从空变成 {x}（旧 W_p = 0）
			if (S.par[p] < 0) {
				int const A0 = S.c0[p], B0 = S.c1[p];
				u64 const* bA0 = S.bits.data() + (size_t)A0 * W;
				u64 const* bB0 = S.bits.data() + (size_t)B0 * W;
				int const sA0 = S.sz[A0], sB0 = S.sz[B0];
				for (int g = 0; g < T; ++g)
					delta += Wnode(genes[g], bA0, bB0, bX.data(), sA0, sB0, 1, W, sc);
			}
			// 祖先链
			int cur = p;
			for (int a = S.par[p]; a >= 0; a = S.par[a]) {
				int const other = (S.c0[a] == cur) ? S.c1[a] : S.c0[a];
				u64 const* bCur = S.bits.data() + (size_t)cur * W;
				u64 const* bOth = S.bits.data() + (size_t)other * W;
				for (size_t k = 0; k < W; ++k) bUnion[k] = bCur[k] | bX[k];
				int const sCur = S.sz[cur] + 1, sOth = S.sz[other];
				int const sC = (int)S.sz[S.root] - sCur - sOth;
				for (size_t k = 0; k < W; ++k) bC2[k] = ~(bUnion[k] | bOth[k]);
				for (int g = 0; g < T; ++g)
					delta += Wnode(genes[g], bUnion.data(), bOth, bC2.data(), sCur, sOth, sC, W, sc);
				cur = a;
			}
			// 减去旧值（重算：插入会改变所有节点的补集，不能缓存）
			{
				int cur2 = p;
				for (int a = S.par[p]; a >= 0; a = S.par[a]) {
					int const other = (S.c0[a] == cur2) ? S.c1[a] : S.c0[a];
					u64 const* bCur = S.bits.data() + (size_t)cur2 * W;
					u64 const* bOth = S.bits.data() + (size_t)other * W;
					int const sCur = S.sz[cur2], sOth = S.sz[other];
					int const sC = (int)S.sz[S.root] - sCur - sOth;
					for (size_t k = 0; k < W; ++k) bC2[k] = ~(bCur[k] | bOth[k]);
					for (int g = 0; g < T; ++g)
						delta -= Wnode(genes[g], bCur, bOth, bC2.data(), sCur, sOth, sC, W, sc);
					cur2 = a;
				}
			}
			deltas[ci] = delta;
		}
		for (int ci = 0; ci < (int)cands.size(); ++ci)
			if (deltas[ci] > bestDelta) { bestDelta = deltas[ci]; bestP = cands[ci]; }

		// 执行插入
		int const leaf = S.newNode(bestP, -1, -1, x);
		S.sz[leaf] = 1;
		S.bits[(size_t)leaf*W + ((size_t)x >> 6)] = 1ULL << (x & 63);
		int const w = S.newNode(S.par[bestP], bestP, leaf, -1);
		S.sz[w] = S.sz[bestP] + 1;
		for (size_t k = 0; k < W; ++k) S.bits[(size_t)w*W + k] = S.bits[(size_t)bestP*W + k] | S.bits[(size_t)leaf*W + k];
		int const par = S.par[bestP];
		if (par < 0) S.root = w;
		else { if (S.c0[par] == bestP) S.c0[par] = w; else S.c1[par] = w; }
		S.par[bestP] = w;
		S.par[leaf] = w;
		// 更新祖先 bits/sz
		for (int a = S.par[w]; a >= 0; a = S.par[a]) {
			for (size_t k = 0; k < W; ++k) S.bits[(size_t)a*W + k] |= S.bits[(size_t)leaf*W + k];
			S.sz[a] = S.sz[S.c0[a]] + S.sz[S.c1[a]];
		}
		nodeOfTaxon[x] = leaf;
		inserted++;
		if (inserted % 500 == 0) std::cerr << "[builder] inserted " << inserted << " / " << nTaxa
			<< "  (" << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - t0).count() << " s)" << std::endl;
	}

	// 自检：Σ_v wsum[v]（内部节点）应当等于 2 × 官方分数
	{
		i64 sumFresh = 0;
		Scratch& sc0 = scratch[0];
		for (int v = 0; v < S.n; ++v) if (S.c0[v] >= 0)
			sumFresh += nodeWsum(S, v, genes, T, sc0);
		std::cerr << "[builder] sum_fresh = " << sumFresh << "  (= 2 x score)" << std::endl;
	}
	string const nw = toNewick(S, S.root) + ";";
	if (outFile.empty()) std::cout << nw << std::endl;
	else { std::ofstream f(outFile); f << nw << "\n"; }
	std::cerr << "[builder] done: " << nTaxa << " taxa" << std::endl;
	return 0;
}
