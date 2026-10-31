// builder: 用"打分插入"构造初始物种树（无根四叶树口径）。
//
// 思路：逐个插入类群；一次计算所有旧内部节点的 dA/dB/dC，树上前缀给出各候选增量。
// 前 --budget 个类群用全候选（树小、便宜），其余只在"伙伴锚点"附近的候选里挑。
// 输出一棵二叉 Newick，交给 solver 做 NNI。
//
// 编译: g++ -std=c++20 -O3 -march=native -fopenmp -I ../aster2/src builder.cpp -o builder

#include "newick.hpp"
#include <algorithm>
#include <chrono>
#include <climits>
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
	vector<u64> qA, qB, qC; // 批量插入计算的分区临时 bitset
	vector<i64> c1, a1, m1, G1, H1, c2, a2, m2, G2, H2;
	void resize(int n, size_t W) {
		c1.resize(n); a1.resize(n); m1.resize(n); G1.resize(n); H1.resize(n);
		c2.resize(n); a2.resize(n); m2.resize(n); G2.resize(n); H2.resize(n);
		t1.assign(W, 0); t2.assign(W, 0);
		qA.assign(W, 0); qB.assign(W, 0); qC.assign(W, 0);
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
	if (sA >= 2 && sB + sC >= 2) { for (size_t w = 0; w < W; ++w) t[w] = B[w] | C[w]; gA = (i64)evalF(g, A, t, sc); }
	if (sB >= 2 && sA + sC >= 2) { for (size_t w = 0; w < W; ++w) t[w] = A[w] | C[w]; gB = (i64)evalF(g, B, t, sc); }
	if (sC >= 2 && sA + sB >= 2) { for (size_t w = 0; w < W; ++w) t[w] = A[w] | B[w]; gC = (i64)evalF(g, C, t, sc); }
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

// 计算节点 v 的 W（对给定的基因子集求和；用于插入后刷新缓存）。
// activeMask 是当前部分树中的类群全集；不能把尚未插入的类群计入 C。
i64 nodeWsum(DynTree const& S, int v, vector<GeneTree> const& genes, int T,
	vector<u64> const& activeMask, int activeCount, Scratch& sc) {
	int const a = S.c0[v], b = S.c1[v];
	if (a < 0) return 0;
	u64 const* bA = S.bits.data() + (size_t)a * S.W;
	u64 const* bB = S.bits.data() + (size_t)b * S.W;
	int const sC = activeCount - S.sz[a] - S.sz[b];
	for (size_t w = 0; w < S.W; ++w)
		sc.t1[w] = activeMask[w] & ~(bA[w] | bB[w]);
	i64 sum = 0;
	for (int g = 0; g < T; ++g)
		sum += Wnode(genes[g], bA, bB, sc.t1.data(), S.sz[a], S.sz[b], sC, S.W, sc);
	return sum;
}

inline void setComplement(u64* dst, u64 const* universe, u64 const* A,
	u64 const* B, size_t W) {
	for (size_t w = 0; w < W; ++w) dst[w] = universe[w] & ~(A[w] | B[w]);
}

// 精确计算把 x 插到 p 上方时 ΣW 的变化。除插入路径上簇发生扩张外，
// 其余旧节点的 C 分区也会增加 x，因此所有旧内部节点都必须重新计分。
i64 candidateDelta(DynTree const& S, int p, u64 const* bX,
	vector<u64> const& activePlusX, vector<GeneTree> const& genes, int T, Scratch& sc) {
	size_t const W = S.W;
	int const afterCount = S.sz[S.root] + 1;
	u64 const* const bP = S.bits.data() + (size_t)p * W;
	vector<int> pathChild(S.n, -1);
	int cur = p;
	for (int a = S.par[p]; a >= 0; a = S.par[a]) {
		pathChild[a] = cur;
		cur = a;
	}

	vector<u64> bA(W), bB(W), cNew(W);
	i64 delta = 0;
	for (int v = 0; v < S.n; ++v) {
		if (S.c0[v] < 0) continue;
		int const left = S.c0[v], right = S.c1[v];
		u64 const* const oldA = S.bits.data() + (size_t)left * W;
		u64 const* const oldB = S.bits.data() + (size_t)right * W;
		bool const expandLeft = pathChild[v] == left;
		bool const expandRight = pathChild[v] == right;
		for (size_t w = 0; w < W; ++w) {
			bA[w] = oldA[w] | (expandLeft ? bX[w] : 0);
			bB[w] = oldB[w] | (expandRight ? bX[w] : 0);
		}
		setComplement(cNew.data(), activePlusX.data(), bA.data(), bB.data(), W);
		int const sizeA = S.sz[left] + (expandLeft ? 1 : 0);
		int const sizeB = S.sz[right] + (expandRight ? 1 : 0);
		i64 newW = 0;
		for (int g = 0; g < T; ++g)
			newW += Wnode(genes[g], bA.data(), bB.data(), cNew.data(), sizeA, sizeB,
				afterCount - sizeA - sizeB, W, sc);
		delta += newW - S.wsum[v];
	}

	// 插入新增的内部节点 (p, x)。
	setComplement(cNew.data(), activePlusX.data(), bP, bX, W);
	for (int g = 0; g < T; ++g)
		delta += Wnode(genes[g], bP, bX, cNew.data(), S.sz[p], 1,
			afterCount - S.sz[p] - 1, W, sc);
	return delta;
}

struct BatchStats {
	vector<i64> dA, dB, dC;       // x 加入旧节点的左、右、外部分区增量
	vector<i64> newNode;          // 候选位置新增内部节点的 W 增量
	vector<i64> prefix;           // 从根到节点的方向修正前缀
	i64 base = 0;                 // Σ_v dC[v]
};

// 一次为所有旧节点计算 dA/dB/dC，并为所有候选位置计算新增节点贡献。
// 旧 candidateDelta 保留作 oracle；本函数只改变评估顺序，不改变 W 定义。
void computeBatchStats(DynTree const& S, u64 const* bX,
	vector<u64> const& activePlusX, vector<int> const& candidates,
	vector<GeneTree> const& genes, int T,
	vector<Scratch>& scratch, int threads, BatchStats& out) {
	size_t const W = S.W;
	int const afterCount = S.sz[S.root] + 1;
	out.dA.assign(S.n, 0);
	out.dB.assign(S.n, 0);
	out.dC.assign(S.n, 0);
	out.newNode.assign(S.n, 0);
	out.prefix.assign(S.n, 0);
	out.base = 0;

	#pragma omp parallel for schedule(dynamic) num_threads(threads)
	for (int v = 0; v < S.n; ++v) {
		if (S.c0[v] < 0) continue;
		Scratch& sc = scratch[omp_get_thread_num()];
		int const left = S.c0[v], right = S.c1[v];
		u64 const* const oldA = S.bits.data() + (size_t)left * W;
		u64 const* const oldB = S.bits.data() + (size_t)right * W;
		i64 const oldW = S.wsum[v];

		// x 加入 A：只扩展 A，C 取 activePlusX 的补集。
		for (size_t w = 0; w < W; ++w) {
			sc.qA[w] = oldA[w] | bX[w];
			sc.qB[w] = oldB[w];
			sc.qC[w] = activePlusX[w] & ~(sc.qA[w] | sc.qB[w]);
		}
		i64 newA = 0;
		for (int g = 0; g < T; ++g)
			newA += Wnode(genes[g], sc.qA.data(), sc.qB.data(), sc.qC.data(),
				S.sz[left] + 1, S.sz[right], afterCount - S.sz[left] - S.sz[right] - 1, W, sc);
		out.dA[v] = newA - oldW;

		// x 加入 B。
		for (size_t w = 0; w < W; ++w) {
			sc.qA[w] = oldA[w];
			sc.qB[w] = oldB[w] | bX[w];
			sc.qC[w] = activePlusX[w] & ~(sc.qA[w] | sc.qB[w]);
		}
		i64 newB = 0;
		for (int g = 0; g < T; ++g)
			newB += Wnode(genes[g], sc.qA.data(), sc.qB.data(), sc.qC.data(),
				S.sz[left], S.sz[right] + 1, afterCount - S.sz[left] - S.sz[right] - 1, W, sc);
		out.dB[v] = newB - oldW;

		// x 加入 C；A/B 不变，C 的大小增加一。
		for (size_t w = 0; w < W; ++w) {
			sc.qA[w] = oldA[w];
			sc.qB[w] = oldB[w];
			sc.qC[w] = activePlusX[w] & ~(sc.qA[w] | sc.qB[w]);
		}
		i64 newC = 0;
		for (int g = 0; g < T; ++g)
			newC += Wnode(genes[g], sc.qA.data(), sc.qB.data(), sc.qC.data(),
				S.sz[left], S.sz[right], afterCount - S.sz[left] - S.sz[right], W, sc);
		out.dC[v] = newC - oldW;
	}

	for (int v = 0; v < S.n; ++v)
		if (S.c0[v] >= 0) out.base += out.dC[v];

	// 节点编号按插入顺序增长，不保证是拓扑序；显式从 root 传播前缀。
	vector<int> st{ S.root };
	while (!st.empty()) {
		int const v = st.back();
		st.pop_back();
		if (S.c0[v] < 0) continue;
		int const left = S.c0[v], right = S.c1[v];
		out.prefix[left] = out.prefix[v] + out.dA[v] - out.dC[v];
		out.prefix[right] = out.prefix[v] + out.dB[v] - out.dC[v];
		st.push_back(left);
		st.push_back(right);
	}

	// 新增节点 (p,x) 的 W 与旧节点方向无关，可按候选并行计算。
	#pragma omp parallel for schedule(dynamic) num_threads(threads)
	for (int ci = 0; ci < (int)candidates.size(); ++ci) {
		int const p = candidates[ci];
		Scratch& sc = scratch[omp_get_thread_num()];
		u64 const* const bP = S.bits.data() + (size_t)p * W;
		for (size_t w = 0; w < W; ++w)
			sc.qC[w] = activePlusX[w] & ~(bP[w] | bX[w]);
		for (int g = 0; g < T; ++g)
			out.newNode[p] += Wnode(genes[g], bP, bX, sc.qC.data(),
				S.sz[p], 1, afterCount - S.sz[p] - 1, W, sc);
	}
}

void insertFirstTaxon(DynTree& S, int x) {
	S.taxon.push_back(x); S.par.push_back(-1); S.c0.push_back(-1); S.c1.push_back(-1);
	S.sz.push_back(1); S.wsum.push_back(0); S.bits.resize(S.W, 0);
	S.bits[(size_t)x >> 6] |= 1ULL << (x & 63);
	S.n = 1; S.root = 0;
}

void insertSecondTaxon(DynTree& S, int x) {
	int const leaf = S.newNode(S.root, -1, -1, x);
	S.bits[(size_t)leaf * S.W + ((size_t)x >> 6)] |= 1ULL << (x & 63);
	S.sz[leaf] = 1;
	int const newRoot = S.newNode(-1, S.root, leaf, -1);
	S.par[S.root] = newRoot; S.par[leaf] = newRoot;
	for (size_t w = 0; w < S.W; ++w)
		S.bits[(size_t)newRoot * S.W + w] = S.bits[(size_t)S.c0[newRoot] * S.W + w]
			| S.bits[(size_t)S.c1[newRoot] * S.W + w];
	S.sz[newRoot] = S.sz[S.c0[newRoot]] + S.sz[S.c1[newRoot]];
	S.root = newRoot;
}

int insertTaxonAbove(DynTree& S, int p, int x) {
	int const leaf = S.newNode(p, -1, -1, x);
	S.sz[leaf] = 1;
	S.bits[(size_t)leaf * S.W + ((size_t)x >> 6)] = 1ULL << (x & 63);
	int const par = S.par[p];
	int const w = S.newNode(par, p, leaf, -1);
	S.sz[w] = S.sz[p] + 1;
	for (size_t k = 0; k < S.W; ++k)
		S.bits[(size_t)w * S.W + k] = S.bits[(size_t)p * S.W + k]
			| S.bits[(size_t)leaf * S.W + k];
	if (par < 0) S.root = w;
	else if (S.c0[par] == p) S.c0[par] = w;
	else S.c1[par] = w;
	S.par[p] = w;
	S.par[leaf] = w;
	for (int a = S.par[w]; a >= 0; a = S.par[a]) {
		for (size_t k = 0; k < S.W; ++k)
			S.bits[(size_t)a * S.W + k] |= S.bits[(size_t)leaf * S.W + k];
		S.sz[a] = S.sz[S.c0[a]] + S.sz[S.c1[a]];
	}
	return leaf;
}

int commitBatchInsertion(DynTree& S, int p, int x, vector<u64>& activeMask,
	BatchStats const& stats) {
	int const oldN = S.n;
	// 必须在改拓扑前保存 p 到根的原始方向。
	vector<std::pair<int, int>> path;
	for (int child = p, a = S.par[p]; a >= 0; child = a, a = S.par[a]) {
		int const dir = (S.c0[a] == child) ? 0 : 1;
		path.emplace_back(a, dir);
	}

	// 所有旧节点先按 x 落入 C 更新；严格祖先再修正为其所在孩子方向。
	for (int v = 0; v < oldN; ++v)
		if (S.c0[v] >= 0) S.wsum[v] += stats.dC[v];
	for (auto const& [a, dir] : path)
		S.wsum[a] += (dir == 0 ? stats.dA[a] : stats.dB[a]) - stats.dC[a];

	int const leaf = insertTaxonAbove(S, p, x);
	activeMask[(size_t)x >> 6] |= 1ULL << (x & 63);
	int const newInternal = oldN + 1; // insertTaxonAbove 先建 leaf，再建 internal
	S.wsum[leaf] = 0;
	S.wsum[newInternal] = stats.newNode[p];
	return leaf;
}

int commitInsertion(DynTree& S, int p, int x, vector<u64>& activeMask,
	vector<GeneTree> const& genes, int T, Scratch& sc) {
	int const leaf = insertTaxonAbove(S, p, x);
	activeMask[(size_t)x >> 6] |= 1ULL << (x & 63);
	int const activeCount = S.sz[S.root];
	for (int v = 0; v < S.n; ++v)
		if (S.c0[v] >= 0)
			S.wsum[v] = nodeWsum(S, v, genes, T, activeMask, activeCount, sc);
	return leaf;
}

string toNewick(DynTree const& S, int u) {
	if (S.taxon[u] >= 0) return common::taxonName2ID[(size_t)S.taxon[u]];
	return "(" + toNewick(S, S.c0[u]) + "," + toNewick(S, S.c1[u]) + ")";
}

} // namespace

#ifndef BUILDER_SELFTEST
int main(int argc, char** argv) {
	string inFile, outFile, seedOrder;
	int threads = 8, budget = 300, cand = 10, Tgenes = 8, orderMode = 0;
	bool nosearch = false;
	long seed = -1;
	for (int i = 1; i < argc; ++i) {
		string s = argv[i];
		auto nx = [&]() -> string { if (i + 1 < argc) return argv[++i]; std::exit(2); };
		if (s == "-i") inFile = nx();
		else if (s == "-o") outFile = nx();
		else if (s == "-t") threads = std::stoi(nx());
		else if (s == "--budget") budget = std::stoi(nx());
		else if (s == "--cand") cand = std::stoi(nx());
		else if (s == "--genes") Tgenes = std::stoi(nx());
		else if (s == "--nosearch") nosearch = true;
		else if (s == "--seed") seed = std::stol(nx());
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
	if (seed >= 0) {
		std::mt19937_64 rng((std::uint64_t)seed);
		std::shuffle(order.begin(), order.end(), rng);
	}

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
	vector<u64> activeMask(W, 0);
	vector<Scratch> scratch(threads);
	int maxN = 1; for (auto const& g : genes) maxN = std::max(maxN, g.n);
	for (auto& sc : scratch) sc.resize(maxN, W);

	vector<int> nodeOfTaxon(nTaxa, -1);   // 已在树中的类群 -> 叶子节点 id
	int inserted = 0;
	auto t0 = std::chrono::steady_clock::now();

	// 以 order 逐个插入
	for (int idx = 0; idx < (int)nTaxa; ++idx) {
		int const x = order[idx];
		if (inserted == 0) {
			insertFirstTaxon(S, x);
			activeMask[(size_t)x >> 6] |= 1ULL << (x & 63);
			nodeOfTaxon[x] = 0;
			inserted = 1;
			continue;
		}
		if (inserted == 1) {
			int const leaf = S.n;
			insertSecondTaxon(S, x);
			activeMask[(size_t)x >> 6] |= 1ULL << (x & 63);
			nodeOfTaxon[x] = leaf;
			inserted = 2;
			continue;
		}

		// 候选位置集合
		vector<int> cands;
		if (nosearch) {
			cands.push_back(S.root);          // 直接挂到根：O(1) 插入，复现"顺序毛毛虫"
		} else if (inserted < budget) {
			for (int v = 0; v < S.n; ++v) cands.push_back(v);
		} else {
			int anchor = anchorTaxon[x];
			int v = (anchor >= 0) ? nodeOfTaxon[anchor] : S.root;
			if (v < 0) v = S.root;
			for (int step = 0; step < cand && v >= 0; ++step) { cands.push_back(v); v = S.par[v]; }
			// 锚点的兄弟孩子
			// （略：路径已覆盖大部分情况）
		}

		// 候选来源扩展后仍保持稳定的首选顺序，并避免重复位置并行写入。
		{
			vector<char> seen(S.n, 0);
			vector<int> unique;
			unique.reserve(cands.size());
			for (int p : cands) if (!seen[p]) {
				seen[p] = 1;
				unique.push_back(p);
			}
			cands.swap(unique);
		}

		// 评估：先一次计算所有旧节点的三个方向增量，再 O(1) 得到候选总增量。
		vector<u64> bX(W, 0), activePlusX = activeMask;
		bX[(size_t)x >> 6] |= 1ULL << (x & 63);
		for (size_t k = 0; k < W; ++k) activePlusX[k] |= bX[k];
		BatchStats stats;
		computeBatchStats(S, bX.data(), activePlusX, cands, genes, T, scratch, threads, stats);
		i64 bestDelta = LLONG_MIN; int bestP = cands.empty() ? S.root : cands[0];
		vector<i64> deltas(cands.size(), LLONG_MIN);
		for (size_t ci = 0; ci < cands.size(); ++ci) {
			int const p = cands[ci];
			deltas[ci] = stats.base + stats.prefix[p] + stats.newNode[p];
			if (deltas[ci] > bestDelta) { bestDelta = deltas[ci]; bestP = p; }
		}

		// 执行插入：旧节点缓存只按方向增量更新，不再全量重建。
		int const leaf = commitBatchInsertion(S, bestP, x, activeMask, stats);
		nodeOfTaxon[x] = leaf;
		++inserted;
		if (inserted % 500 == 0) std::cerr << "[builder] inserted " << inserted << " / " << nTaxa
			<< "  (" << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - t0).count() << " s)" << std::endl;
	}

	// 自检：Σ_v wsum[v]（内部节点）应当等于 2 × 官方分数
	{
		i64 sumFresh = 0;
		Scratch& sc0 = scratch[0];
		for (int v = 0; v < S.n; ++v) if (S.c0[v] >= 0)
			sumFresh += nodeWsum(S, v, genes, T, activeMask, inserted, sc0);
		i64 sumCached = 0;
		for (int v = 0; v < S.n; ++v) if (S.c0[v] >= 0) sumCached += S.wsum[v];
		if (sumCached != sumFresh) {
			std::cerr << "[builder] ERROR: cached sum " << sumCached
				<< " != fresh sum " << sumFresh << std::endl;
			return 1;
		}
		std::cerr << "[builder] sum_fresh = " << sumFresh << " (= 2 x score)" << std::endl;
	}
	string const nw = toNewick(S, S.root) + ";";
	if (outFile.empty()) std::cout << nw << std::endl;
	else { std::ofstream f(outFile); f << nw << "\n"; }
	std::cerr << "[builder] done: " << nTaxa << " taxa" << std::endl;
	return 0;
}
#endif
