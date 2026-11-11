// nni_incr.cpp — 增量 13 状态 DP 驱动的 NNI 局部搜索（替代 evalF 六项分解）。
//
// 原理：
//   分数 = Σ_v W_v / 2，W_v 只依赖节点 v 自己的三分划 (A,B,C)；
//   NNI 只改 u、v 两个节点的孩子簇 ⇒ 只需重算这两个 W。
//   每个基因树先 Huffman resolve 成二叉（合成边 C），维护 V1/V2 两套 13 状态；
//   改一个叶子的颜色只重算它到根的祖先链 ⇒ O(深度)。
//
// 关键实现点（勿改坏）：
//   1. 物种节点按前序枚举，染色状态跨节点滚动复用（代价 = |旧簇| + |新簇|）；
//   2. 每个候选变体的染色顺序固定为 S0(当前 v) → S1(v') → S2(u') → 回滚 S0；
//   3. 物种节点被改动时只有 v 的叶集合变化（u 及所有祖先的簇集合不变），
//      所以只需重建 taxa[v]。
//
// 构建: g++ -std=c++20 -O3 -march=native -fopenmp nni_incr.cpp -o nni_incr
//
// 用法:
//   ./nni_incr -s sp.tree -i genes.trees -o out.tree [--sweeps N] [-t P] [--genes N]
//   ./nni_incr -s sp.tree -i genes.trees --score      # 只算初始分
//   ./nni_incr -s sp.tree -i genes.trees --check 10   # 增量 vs 从零重算 对拍
#include <array>
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <queue>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>
#include <omp.h>

using ll = long long;
using St = std::array<ll, 13>;

static double now() {
	return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// ------------------------------------------------------------------ 解析
struct Raw { std::vector<int> ch; std::string lab, name; };

static std::vector<Raw> parseRaw(const std::string& s0, int& root) {
	std::string s;
	s.reserve(s0.size());
	for (char c : s0) if (c != '\r' && c != '\n') s.push_back(c);
	std::vector<Raw> P; std::vector<int> st; root = -1;
	int n = (int)s.size(), i = 0;
	auto isD = [](char c) { return c == ',' || c == '(' || c == ')' || c == ';'; };
	while (i < n) {
		char c = s[i];
		if (c == '(') {
			int id = (int)P.size(); P.push_back({});
			if (!st.empty()) P[st.back()].ch.push_back(id); else root = id;
			st.push_back(id); ++i;
		}
		else if (c == ',') ++i;
		else if (c == ')') {
			++i; int j = i; while (j < n && !isD(s[j])) ++j;
			if (!st.empty()) { P[st.back()].lab = s.substr(i, j - i); st.pop_back(); }
			i = j;
		}
		else if (c == ';') ++i;
		else {
			int j = i; while (j < n && !isD(s[j])) ++j;
			int id = (int)P.size(); P.push_back({}); P[id].name = s.substr(i, j - i);
			if (!st.empty()) P[st.back()].ch.push_back(id); else root = id;
			i = j;
		}
	}
	return P;
}

// ------------------------------------------------- 基因树（Huffman 二叉化 + C 边）
struct Bin {
	std::vector<int> lc, rc, par, sz, tid;
	std::vector<unsigned char> code;   // 3=叶 0=C(合成) 1=S 2=D
	std::vector<int> leaves, order;    // order = 后序
	std::vector<int> leafOfTaxon;      // size nTaxa，该基因中该类的叶节点 id 或 -1
	std::vector<unsigned char> col;    // 叶子的当前颜色
	std::vector<St> V1, V2;
	int root = 0, n = 0;
};

static Bin toBinary(const std::vector<Raw>& P, int root,
	const std::unordered_map<std::string, int>& name2id, int nTaxa) {
	Bin b;
	std::function<int(int)> conv = [&](int u) -> int {
		int id = b.n++;
		b.lc.push_back(-1); b.rc.push_back(-1); b.par.push_back(-1); b.sz.push_back(0);
		b.code.push_back(P[u].ch.empty() ? 3 : (P[u].lab == "C" ? 0 : (P[u].lab == "D" ? 2 : 1)));
		b.tid.push_back(-1);
		if (P[u].ch.empty()) {
			auto it = name2id.find(P[u].name);
			b.tid[id] = (it == name2id.end()) ? -1 : it->second;
			b.leaves.push_back(id); b.sz[id] = 1;
			return id;
		}
		std::vector<int> ids;
		for (int c : P[u].ch) ids.push_back(conv(c));
		if (ids.size() == 1) {                     // 单叉：折叠到唯一孩子
			int c = ids[0];
			b.lc[id] = b.lc[c]; b.rc[id] = b.rc[c];
			b.code[id] = b.code[c]; b.tid[id] = b.tid[c]; b.sz[id] = b.sz[c];
			if (b.lc[id] >= 0) { b.par[b.lc[id]] = id; b.par[b.rc[id]] = id; }
			else b.leaves.back() = id;
			return id;
		}
		if (ids.size() == 2) { b.lc[id] = ids[0]; b.rc[id] = ids[1]; }
		else {                                     // 多叉：按子树大小 Huffman 合并
			std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>,
				std::greater<std::pair<int, int>>> pq;
			for (int c : ids) pq.push({ b.sz[c], c });
			while (pq.size() > 2) {
				auto [s1, x] = pq.top(); pq.pop();
				auto [s2, y] = pq.top(); pq.pop();
				int nid = b.n++;
				b.lc.push_back(x); b.rc.push_back(y); b.par.push_back(-1);
				b.sz.push_back(s1 + s2); b.code.push_back(0); b.tid.push_back(-1);
				b.par[x] = nid; b.par[y] = nid;
				pq.push({ s1 + s2, nid });
			}
			auto [s1, x] = pq.top(); pq.pop();
			auto [s2, y] = pq.top(); pq.pop();
			b.lc[id] = x; b.rc[id] = y;
		}
		b.par[b.lc[id]] = id; b.par[b.rc[id]] = id;
		b.sz[id] = b.sz[b.lc[id]] + b.sz[b.rc[id]];
		return id;
	};
	b.root = conv(root);
	b.n = (int)b.lc.size();
	std::vector<int> st{ b.root };
	while (!st.empty()) {
		int t = st.back(); st.pop_back();
		b.order.push_back(t);
		if (b.lc[t] >= 0) { st.push_back(b.lc[t]); st.push_back(b.rc[t]); }
	}
	std::reverse(b.order.begin(), b.order.end());
	b.leafOfTaxon.assign(nTaxa, -1);
	b.col.assign(b.n, 2);
	for (int l : b.leaves) if (b.tid[l] >= 0) b.leafOfTaxon[b.tid[l]] = l;
	b.V1.assign(b.n, St{}); b.V2.assign(b.n, St{});
	for (auto& s : b.V1) s.fill(0);
	for (auto& s : b.V2) s.fill(0);
	return b;
}

// ------------------------------------------------------------ 13 状态步进（勿改）
static inline void step1(St& A, const St& U, const St& V, int nu, int nv) {
	auto [Ru, Bu, Gu, RRu, BBu, GUu, RBu, RGu, BGu, coRu, coBu, coGu, Qu] = U;
	auto [Rv, Bv, Gv, RRv, BBv, GUv, RBv, RGv, BGv, coRv, coBv, coGv, Qv] = V;
	A[0] = nu * Ru + nv * Rv; A[1] = nu * Bu + nv * Bv; A[2] = nu * Gu + nv * Gv;
	A[3] = RRu + RRv + (ll)nu * nv * Ru * Rv;
	A[4] = BBu + BBv + (ll)nu * nv * Bu * Bv;
	A[5] = GUu + GUv + (ll)nu * nv * Gu * Gv;
	A[6] = RBu + RBv + (ll)nu * nv * (Ru * Bv + Bu * Rv);
	A[7] = RGu + RGv + (ll)nu * nv * (Ru * Gv + Gu * Rv);
	A[8] = BGu + BGv + (ll)nu * nv * (Bu * Gv + Gu * Bv);
	A[9] = nu * coRu + nv * coRv + (BBu * nv * Gv + BBv * nu * Gu + GUu * nv * Bv + GUv * nu * Bu + BGu * nv * Rv + BGv * nu * Ru);
	A[10] = nu * coBu + nv * coBv + (RRu * nv * Gv + RRv * nu * Gu + GUu * nv * Rv + GUv * nu * Ru + RGu * nv * Bv + RGv * nu * Bu);
	A[11] = nu * coGu + nv * coGv + (BBu * nv * Rv + BBv * nu * Ru + RRu * nv * Bv + RRv * nu * Bu + RBu * nv * Gv + RBv * nu * Gu);
	A[12] = Qu + Qv + (ll)nu * nv * (coRu * Rv + coRv * Ru + coBu * Bv + coBv * Bu + coGu * Gv + coGv * Gu)
		+ (RRu * BGv + RRv * BGu + BBu * RGv + BBv * RGu + GUu * RBv + GUv * RBu);
}

static inline void step2(St& A, const St& U, const St& V, int cu, int cv, int nu, int nv) {
	auto [Ru, Bu, Gu, RRu, BBu, GUu, RBu, RGu, BGu, coRu, coBu, coGu, Qu] = U;
	auto [Rv, Bv, Gv, RRv, BBv, GUv, RBv, RGv, BGv, coRv, coBv, coGv, Qv] = V;
	A[0] = nu * Ru + nv * Rv; A[1] = nu * Bu + nv * Bv; A[2] = nu * Gu + nv * Gv;
	A[3] = cu * RRu + cv * RRv + (ll)nu * nv * Ru * Rv;
	A[4] = cu * BBu + cv * BBv + (ll)nu * nv * Bu * Bv;
	A[5] = cu * GUu + cv * GUv + (ll)nu * nv * Gu * Gv;
	A[6] = cu * RBu + cv * RBv + (ll)nu * nv * (Ru * Bv + Bu * Rv);
	A[7] = cu * RGu + cv * RGv + (ll)nu * nv * (Ru * Gv + Gu * Rv);
	A[8] = cu * BGu + cv * BGv + (ll)nu * nv * (Bu * Gv + Gu * Bv);
	A[9] = nu * coRu + nv * coRv + (ll)cu * nv * (BBu * Gv + GUu * Bv + BGu * Rv) + (ll)cv * nu * (BBv * Gu + GUv * Bu + BGv * Ru);
	A[10] = nu * coBu + nv * coBv + (ll)cu * nv * (RRu * Gv + GUu * Rv + RGu * Bv) + (ll)cv * nu * (RRv * Gu + GUv * Ru + RGv * Bu);
	A[11] = nu * coGu + nv * coGv + (ll)cu * nv * (RRu * Bv + BBu * Rv + RBu * Gv) + (ll)cv * nu * (RRv * Bu + BBv * Ru + RBv * Gu);
	A[12] = Qu + Qv + (ll)nu * nv * (coRu * Rv + coRv * Ru + coBu * Bv + coBv * Bu + coGu * Gv + coGv * Gu)
		+ (ll)cu * cv * (RRu * BGv + RRv * BGu + BBu * RGv + BBv * RGu + GUu * RBv + GUv * RBu);
}

static inline void recomputeNode(Bin& b, int a, int l, int r) {
	int nu = (b.code[l] == 2) ? 0 : 1, nv = (b.code[r] == 2) ? 0 : 1;
	int cu = (b.code[l] == 0) ? 1 : 0, cv = (b.code[r] == 0) ? 1 : 0;
	step1(b.V1[a], b.V1[l], b.V1[r], nu, nv);
	step2(b.V2[a], b.V2[l], b.V2[r], cu, cv, nu, nv);
}

static inline void recolorLeaf(Bin& b, int leaf, int c) {
	if (b.col[leaf] == (unsigned char)c) return;
	b.col[leaf] = (unsigned char)c;
	b.V1[leaf].fill(0); b.V2[leaf].fill(0);
	b.V1[leaf][c] = 1; b.V2[leaf][c] = 1;
	for (int a = b.par[leaf]; a >= 0; a = b.par[a]) recomputeNode(b, a, b.lc[a], b.rc[a]);
}

static inline ll readW(const Bin& b) { return b.V1[b.root][12] - b.V2[b.root][12]; }

// ------------------------------------------------------------ 物种树
struct Sp {
	std::vector<int> par, c0, c1;
	std::vector<std::vector<int>> taxa;   // 每个节点子树的类群（只保证集合正确）
	std::vector<int> internal;            // 前序
	std::vector<std::string> id2name;
	int n = 0, root = 0, nTaxa = 0;
};

static inline void paintCluster(Bin& b, const std::vector<int>& ts, int c) {
	for (int t : ts) { int l = b.leafOfTaxon[t]; if (l >= 0) recolorLeaf(b, l, c); }
}

// 旧三分划 -> 新三分划（oldA<0 表示"无旧三分划"）
static inline void setPartition(Bin& b, const Sp& sp, int oldA, int oldB, int newA, int newB) {
	if (oldA >= 0) { paintCluster(b, sp.taxa[oldA], 2); paintCluster(b, sp.taxa[oldB], 2); }
	paintCluster(b, sp.taxa[newA], 0);
	paintCluster(b, sp.taxa[newB], 1);
}

// 从任意染色恢复到 (newA,newB) 三分划（校验模式用）
static inline void resetPartition(Bin& b, const Sp& sp, int newA, int newB) {
	for (int l : b.leaves) recolorLeaf(b, l, 2);
	paintCluster(b, sp.taxa[newA], 0);
	paintCluster(b, sp.taxa[newB], 1);
}

// 从零重算某三分划的 W（不用增量，校验用）
static ll fullW(Bin& b, const std::vector<int>& X, const std::vector<int>& Y) {
	for (int l : b.leaves) { b.col[l] = 2; b.V1[l].fill(0); b.V2[l].fill(0); b.V1[l][2] = 1; b.V2[l][2] = 1; }
	for (int t : X) { int l = b.leafOfTaxon[t]; if (l >= 0) { b.col[l] = 0; b.V1[l].fill(0); b.V2[l].fill(0); b.V1[l][0] = 1; b.V2[l][0] = 1; } }
	for (int t : Y) { int l = b.leafOfTaxon[t]; if (l >= 0) { b.col[l] = 1; b.V1[l].fill(0); b.V2[l].fill(0); b.V1[l][1] = 1; b.V2[l][1] = 1; } }
	for (int t : b.order) if (b.lc[t] >= 0) recomputeNode(b, t, b.lc[t], b.rc[t]);
	return readW(b);
}

// ------------------------------------------------------------------ 主程序
int main(int argc, char** argv) {
	std::string spFile, geneFile, outFile;
	int threads = 8, sweeps = 20, maxGenes = 0;
	bool scoreOnly = false, verbose = true;
	bool acceptZero = false, shuffle = false;
	int nCheck = 0, tabuLen = 3, patience = 4, kick = 0, kickEvery = 5;
	unsigned seed = 12345;
	for (int i = 1; i < argc; ++i) {
		std::string s = argv[i];
		auto nx = [&]() -> std::string { if (i + 1 < argc) return argv[++i]; std::exit(2); };
		if (s == "-s") spFile = nx();
		else if (s == "-i") geneFile = nx();
		else if (s == "-o") outFile = nx();
		else if (s == "-t") threads = std::stoi(nx());
		else if (s == "--sweeps") sweeps = std::stoi(nx());
		else if (s == "--genes") maxGenes = std::stoi(nx());
		else if (s == "--score") scoreOnly = true;
		else if (s == "--check") nCheck = std::stoi(nx());
		else if (s == "--zero") acceptZero = true;
		else if (s == "--shuffle") shuffle = true;
		else if (s == "--tabu") tabuLen = std::stoi(nx());
		else if (s == "--patience") patience = std::stoi(nx());
		else if (s == "--seed") seed = (unsigned)std::stoul(nx());
		else if (s == "--kick") kick = std::stoi(nx());
		else if (s == "--kick-every") kickEvery = std::stoi(nx());
		else { std::fprintf(stderr, "unknown arg %s\n", s.c_str()); return 2; }
	}
	if (spFile.empty() || geneFile.empty()) { std::fprintf(stderr, "usage: -s sp -i genes\n"); return 2; }

	std::vector<std::string> spL, gL;
	{
		std::ifstream f(spFile); std::string l;
		while (std::getline(f, l)) if (l.find_first_not_of(" \t\r\n") != std::string::npos) spL.push_back(l);
	}
	{
		std::ifstream f(geneFile); std::string l;
		while (std::getline(f, l)) if (l.find_first_not_of(" \t\r\n") != std::string::npos) gL.push_back(l);
	}
	if (spL.empty() || gL.empty()) { std::fprintf(stderr, "empty input\n"); return 1; }
	if (maxGenes > 0 && (int)gL.size() > maxGenes) gL.resize(maxGenes);

	// --- 物种树
	int spRoot;
	std::vector<Raw> spP = parseRaw(spL[0], spRoot);
	for (auto& nd : spP) if (!nd.ch.empty() && nd.ch.size() != 2) {
		std::fprintf(stderr, "Error: species tree must be binary\n"); return 1;
	}
	Sp sp;
	sp.n = (int)spP.size(); sp.root = spRoot;
	sp.par.assign(sp.n, -1); sp.c0.assign(sp.n, -1); sp.c1.assign(sp.n, -1);
	sp.taxa.assign(sp.n, {});
	std::unordered_map<std::string, int> name2id;
	for (auto& nd : spP) if (nd.ch.empty())
		if (!name2id.count(nd.name)) { name2id[nd.name] = (int)sp.id2name.size(); sp.id2name.push_back(nd.name); }
	sp.nTaxa = (int)sp.id2name.size();
	for (int u = 0; u < sp.n; ++u) {
		if (spP[u].ch.empty()) sp.taxa[u] = { name2id[spP[u].name] };
		else {
			sp.c0[u] = spP[u].ch[0]; sp.c1[u] = spP[u].ch[1];
			sp.par[sp.c0[u]] = u; sp.par[sp.c1[u]] = u;
			sp.internal.push_back(u);
		}
	}
	{
		std::vector<int> ord; std::vector<int> st{ sp.root };
		while (!st.empty()) {
			int t = st.back(); st.pop_back(); ord.push_back(t);
			if (sp.c0[t] >= 0) { st.push_back(sp.c0[t]); st.push_back(sp.c1[t]); }
		}
		std::reverse(ord.begin(), ord.end());
		for (int u : ord) if (sp.c0[u] >= 0) {
			sp.taxa[u] = sp.taxa[sp.c0[u]];
			sp.taxa[u].insert(sp.taxa[u].end(), sp.taxa[sp.c1[u]].begin(), sp.taxa[sp.c1[u]].end());
		}
	}
	{
		// --- 基因树并行二叉化
		const int T = (int)gL.size();
		std::vector<Bin> genes(T);
		double t0 = now();
		#pragma omp parallel for schedule(dynamic, 1) num_threads(threads)
		for (int gi = 0; gi < T; ++gi) {
			int gr;
			std::vector<Raw> gp = parseRaw(gL[gi], gr);
			genes[gi] = toBinary(gp, gr, name2id, sp.nTaxa);
		}
		if (verbose) std::fprintf(stderr, "[incr] genes=%d taxa=%d species_internal=%d  二叉化 %.1fs\n",
			T, sp.nTaxa, (int)sp.internal.size(), now() - t0);

		const int K = (int)sp.internal.size();
		std::vector<int> idx(sp.n, -1);
		for (int k = 0; k < K; ++k) idx[sp.internal[k]] = k;

		// --- 初始评分（前序滚动染色）
		std::vector<ll> cache((size_t)K * T, 0);
		int prevA = -1, prevB = -1;
		double t1 = now();
		for (int k = 0; k < K; ++k) {
			int v = sp.internal[k];
			#pragma omp parallel for schedule(dynamic, 4) num_threads(threads)
			for (int g = 0; g < T; ++g) setPartition(genes[g], sp, prevA, prevB, sp.c0[v], sp.c1[v]);
			prevA = sp.c0[v]; prevB = sp.c1[v];
			#pragma omp parallel for schedule(dynamic, 4) num_threads(threads)
			for (int g = 0; g < T; ++g) cache[(size_t)k * T + g] = readW(genes[g]);
		}
		double tInit = now() - t1;
		__int128 sum = 0;
		for (ll x : cache) sum += x;
		ll total = (ll)(sum / 2);
		if (verbose) std::fprintf(stderr, "[incr] initial score = %lld  (%.2fs)\n", total, tInit);
		if (scoreOnly) { std::printf("%lld\n", total); return 0; }

		// ---------------- 校验模式：增量 vs 从零重算 ----------------
		if (nCheck > 0) {
			int nCmp = 0, nBad = 0; ll maxDiff = 0;
			std::vector<int> uni;
			int step = std::max(1, K / nCheck);
			for (int k = 0; k < K; ++k) {
				int v = sp.internal[k];
				int u = sp.par[v];
				if (u < 0 || sp.c0[u] < 0) continue;
				#pragma omp parallel for schedule(dynamic, 4) num_threads(threads)
				for (int g = 0; g < T; ++g) setPartition(genes[g], sp, prevA, prevB, sp.c0[v], sp.c1[v]);
				prevA = sp.c0[v]; prevB = sp.c1[v];
				if (k % step != 0) continue;
				int a = (sp.c0[u] == v) ? sp.c1[u] : sp.c0[u];
				int c = sp.c0[v], d = sp.c1[v];
				for (int var = 0; var < 2; ++var) {
					int x = (var == 0) ? c : d, y = (var == 0) ? d : c;
					uni = sp.taxa[a];
					uni.insert(uni.end(), sp.taxa[y].begin(), sp.taxa[y].end());
					for (int g = 0; g < std::min(T, 3); ++g) {
						Bin& b = genes[g];
						ll wv, wu;
						if (var == 0) {
							paintCluster(b, sp.taxa[a], 0); paintCluster(b, sp.taxa[c], 2);
							wv = readW(b);
							paintCluster(b, sp.taxa[d], 0); paintCluster(b, sp.taxa[c], 1);
							wu = readW(b);
							paintCluster(b, sp.taxa[a], 2); paintCluster(b, sp.taxa[d], 1); paintCluster(b, sp.taxa[c], 0);
						}
						else {
							paintCluster(b, sp.taxa[a], 0); paintCluster(b, sp.taxa[c], 1); paintCluster(b, sp.taxa[d], 2);
							wv = readW(b);
							paintCluster(b, sp.taxa[c], 0); paintCluster(b, sp.taxa[d], 1);
							wu = readW(b);
							paintCluster(b, sp.taxa[a], 2);
						}
						ll fv = fullW(b, sp.taxa[a], sp.taxa[y]);
						ll fu = fullW(b, uni, sp.taxa[x]);
						resetPartition(b, sp, sp.c0[v], sp.c1[v]);
						++nCmp;
						if (wv != fv || wu != fu) {
							++nBad;
							maxDiff = std::max({ maxDiff, std::llabs(wv - fv), std::llabs(wu - fu) });
							std::fprintf(stderr, "[check] MISMATCH k=%d var=%d g=%d  incr=(%lld,%lld) full=(%lld,%lld)\n",
								k, var, g, wv, wu, fv, fu);
						}
					}
				}
			}
			std::fprintf(stderr, "[check] %d comparisons, %d mismatches, max|diff|=%lld\n", nCmp, nBad, maxDiff);
			return nBad ? 3 : 0;
		}

		// ---------------- NNI 扫描 ----------------
		std::vector<ll> Wv(T), Wu(T);
		std::vector<int> order = sp.internal;
		std::vector<int> tabuUntil(sp.n, 0);
		std::mt19937 rng(seed);
		int noGain = 0;

		// 把染色推进到节点 v 的三分划
		auto prepPartition = [&](int v) {
			#pragma omp parallel for schedule(dynamic, 4) num_threads(threads)
			for (int g = 0; g < T; ++g) setPartition(genes[g], sp, prevA, prevB, sp.c0[v], sp.c1[v]);
			prevA = sp.c0[v]; prevB = sp.c1[v];
		};
		// 评估 (v,var)：写 Wv/Wu，回滚到 v 的三分划，返回 delta
		auto evalVar = [&](int v, int var) -> ll {
			int u = sp.par[v];
			int a = (sp.c0[u] == v) ? sp.c1[u] : sp.c0[u];
			int c = sp.c0[v], d = sp.c1[v];
			int ku = idx[u], kv = idx[v];
			ll delta = 0;
			#pragma omp parallel for schedule(dynamic, 4) num_threads(threads) reduction(+:delta)
			for (int g = 0; g < T; ++g) {
				Bin& b = genes[g];
				ll wv, wu;
				if (var == 0) {
					paintCluster(b, sp.taxa[a], 0); paintCluster(b, sp.taxa[c], 2);
					wv = readW(b);
					paintCluster(b, sp.taxa[d], 0); paintCluster(b, sp.taxa[c], 1);
					wu = readW(b);
					paintCluster(b, sp.taxa[a], 2); paintCluster(b, sp.taxa[d], 1); paintCluster(b, sp.taxa[c], 0);
				}
				else {
					paintCluster(b, sp.taxa[a], 0); paintCluster(b, sp.taxa[c], 1); paintCluster(b, sp.taxa[d], 2);
					wv = readW(b);
					paintCluster(b, sp.taxa[c], 0); paintCluster(b, sp.taxa[d], 1);
					wu = readW(b);
					paintCluster(b, sp.taxa[a], 2);
				}
				Wv[g] = wv; Wu[g] = wu;
				delta += wu + wv - cache[(size_t)ku * T + g] - cache[(size_t)kv * T + g];
			}
			return delta;
		};
		// 应用 (v,var)（调用前必须已 evalVar）
		auto applyVar = [&](int v, int var) {
			int u = sp.par[v];
			int a = (sp.c0[u] == v) ? sp.c1[u] : sp.c0[u];
			int c = sp.c0[v], d = sp.c1[v];
			int x = (var == 0) ? c : d, y = (var == 0) ? d : c;
			int ku = idx[u], kv = idx[v];
			sp.c0[v] = a; sp.c1[v] = y; sp.par[a] = v; sp.par[y] = v;
			sp.c0[u] = v; sp.c1[u] = x; sp.par[v] = u; sp.par[x] = u;
			sp.taxa[v].assign(sp.taxa[a].begin(), sp.taxa[a].end());
			sp.taxa[v].insert(sp.taxa[v].end(), sp.taxa[y].begin(), sp.taxa[y].end());
			#pragma omp parallel for schedule(dynamic, 4) num_threads(threads)
			for (int g = 0; g < T; ++g) {
				cache[(size_t)kv * T + g] = Wv[g];
				cache[(size_t)ku * T + g] = Wu[g];
			}
			// 回滚后染色仍对应"旧 v 的两个孩子"，把 prev 记成它们即可，不必重染
			prevA = c; prevB = d;
		};

		// 最优解快照（ILS 用）
		struct Snap {
			std::vector<int> c0, c1, par;
			std::vector<std::vector<int>> taxa;
			std::vector<ll> cache;
			ll total = 0; bool valid = false;
		} best;
		best.total = total;
		auto saveSnap = [&]() {
			best.c0 = sp.c0; best.c1 = sp.c1; best.par = sp.par;
			best.taxa = sp.taxa; best.cache = cache; best.total = total; best.valid = true;
		};
		auto loadSnap = [&]() {
			sp.c0 = best.c0; sp.c1 = best.c1; sp.par = best.par; sp.taxa = best.taxa;
			cache = best.cache; total = best.total;
			prevA = prevB = -1;
			#pragma omp parallel for schedule(dynamic, 4) num_threads(threads)
			for (int g = 0; g < T; ++g)
				for (int l : genes[g].leaves) recolorLeaf(genes[g], l, 2);
		};

		for (int sweep = 1; sweep <= sweeps; ++sweep) {
			if (kick > 0 && kickEvery > 0 && sweep > 1 && (sweep - 1) % kickEvery == 0) {
				for (int t = 0; t < kick; ++t) {
					int v = order[rng() % order.size()];
					if (sp.par[v] < 0 || sp.c0[sp.par[v]] < 0) continue;
					prepPartition(v);
					int var = (int)(rng() & 1u);
					ll d = evalVar(v, var);
					applyVar(v, var);
					total += d / 2;
				}
				noGain = 0;
			}
			if (shuffle) {
				std::mt19937 rg2(seed + (unsigned)sweep * 7919u);
				std::shuffle(order.begin(), order.end(), rg2);
			}
			ll sweepGain = 0; int posMoves = 0, zeroMoves = 0;
			double ts = now();
			for (int k = 0; k < K; ++k) {
				int v = order[k];
				int u = sp.par[v];
				if (u < 0 || sp.c0[u] < 0) continue;
				prepPartition(v);
				for (int var = 0; var < 2; ++var) {
					ll delta = evalVar(v, var);
					bool acc = (delta > 0) ||
						(acceptZero && delta == 0 && sweep >= tabuUntil[v]);
					if (acc) {
						applyVar(v, var);
						total += delta / 2; sweepGain += delta / 2;
						if (delta == 0) { tabuUntil[v] = sweep + tabuLen; ++zeroMoves; }
						else ++posMoves;
						break;
					}
				}
			}
			if (posMoves > 0) noGain = 0; else ++noGain;
			std::fprintf(stderr, "[incr] sweep %d: +%d =0:%d gain=%lld total=%lld  (%.2fs)\n",
				sweep, posMoves, zeroMoves, sweepGain, total, now() - ts);
			if (total > best.total) { best.total = total; saveSnap(); }
			bool done = (posMoves == 0 && zeroMoves == 0) || (noGain >= patience);
			if (done && kick == 0) break;
		}
		if (best.valid && total < best.total) {
			std::fprintf(stderr, "[incr] ILS: 回退到最优解 %lld（当前 %lld）\n", best.total, total);
			loadSnap();
		}

		// --- 输出
		std::string out; out.reserve(80000);
		std::function<void(int)> rec = [&](int u) {
			if (sp.c0[u] < 0) out += sp.id2name[sp.taxa[u][0]];
			else { out += "("; rec(sp.c0[u]); out += ","; rec(sp.c1[u]); out += ")"; }
		};
		rec(sp.root); out += ";";
		if (outFile.empty()) std::printf("%s\n", out.c_str());
		else { std::ofstream f(outFile); f << out << "\n"; }
		std::fprintf(stderr, "[incr] final score = %lld\n", total);
	}
	return 0;
}
