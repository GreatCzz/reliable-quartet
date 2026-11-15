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
#include "w3.hpp"
#include "sector.hpp"

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

struct Args {
	string species, input, output;
	int threads = 8, rounds = 50, iters = 2000;
	double target = 0.0;   // >0：只从"贡献最低的 target 比例"节点里采样 x
	int screen = 0;        // >0：先用 screen 棵树筛选移动，通过后再用全部树确认
	int localDepth = 0;    // >0：v 优先在 x 原位置附近 localDepth 条边内选（局部 SPR）
	int kick = 0;          // >0：前 kick 次移动无条件接受（扰动重启）
	int allowLeaf = 0;     // 1：允许移动单个叶子（leaf SPR）
	int sectorK = 0;       // >0：每轮先做 K 个类群的 sector 精确重插入
	int sectorGenes = 20;  // sector 排序用的基因数
	string probe;          // "xName:vName" 调试：打印把 x 移到 v 上方的 Δ
	long long seed = 1;
	bool verbose = false;
};

Args parseArgs(int argc, char** argv) {
	Args a;
	for (int i = 1; i < argc; ++i) {
		string s = argv[i];
		auto next = [&]() -> string { if (i + 1 < argc) return argv[++i]; std::exit(2); };
		if (s == "-s") a.species = next();
		else if (s == "-i") a.input = next();
		else if (s == "-o") a.output = next();
		else if (s == "-t") a.threads = std::stoi(next());
		else if (s == "--rounds") a.rounds = std::stoi(next());
		else if (s == "--iters") a.iters = std::stoi(next());
		else if (s == "--seed") a.seed = std::stoll(next());
		else if (s == "--target") a.target = std::stod(next());
		else if (s == "--screen") a.screen = std::stoi(next());
		else if (s == "--local") a.localDepth = std::stoi(next());
		else if (s == "--kick") a.kick = std::stoi(next());
		else if (s == "--leaf") a.allowLeaf = std::stoi(next());
		else if (s == "--sector") a.sectorK = std::stoi(next());
		else if (s == "--sector-genes") a.sectorGenes = std::stoi(next());
		else if (s == "--probe") a.probe = next();
		else if (s == "-v") a.verbose = true;
		else { std::cerr << "unknown arg " << s << std::endl; std::exit(2); }
	}
	return a;
}

} // namespace

int main(int argc, char** argv) {
	Args const args = parseArgs(argc, argv);
	if (args.species.empty() || args.input.empty()) {
		std::cerr << "Error: -s <species tree> and -i <gene trees> are required" << std::endl;
		return 2;
	}
	newick::Newick spTree;
	{
		std::ifstream fin(args.species);
		if (!fin.is_open()) { std::cerr << "Error: cannot open species tree: " << args.species << std::endl; return 2; }
		string line;
		while (std::getline(fin, line)) {
			if (line.find_first_not_of(" \t\r\n") != string::npos) { spTree = newick::Newick(line); break; }
		}
	}
	vector<GeneTree> genes;
	{
		std::ifstream fin(args.input);
		if (!fin.is_open()) { std::cerr << "Error: cannot open gene trees: " << args.input << std::endl; return 2; }
		string line;
		while (std::getline(fin, line)) {
			if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
			newick::Newick nw(line);
			genes.push_back(buildGene(nw.root()));
		}
	}
	int const T = (int)genes.size();
	vector<w3::GeneView> gv;
	for (auto& g : genes) {
		w3::GeneView v;
		v.start = g.start; v.kids = g.kids; v.solid = g.solid; v.taxon = g.taxon; v.root = g.root; v.n = g.n;
		gv.push_back(std::move(v));
	}
	size_t const nTaxa = common::taxonName2ID.nTaxa();
	SpeciesTree S = buildSpecies(spTree.root(), nTaxa);
	int const K = (int)S.internal.size();
	std::cerr << "[spr] nodes=" << S.n << " internal=" << K << " genes=" << T << " taxa=" << nTaxa << std::endl;

	vector<i64> cache((size_t)K * (size_t)T, 0);
	int const P = std::max(1, args.threads);
	vector<Scratch> scratch(P);
	int maxN = 1;
	for (auto const& g : genes) maxN = std::max(maxN, g.n);
	for (auto& s : scratch) { s.resize(maxN); s.t1.assign(S.W, 0); s.t2.assign(S.W, 0); }
	vector<w3::Scratch3> s3(P);
	for (auto& x : s3) x.resize(maxN);

	vector<int> idx(S.n, -1);
	for (int k = 0; k < K; ++k) idx[S.internal[k]] = k;

	// 全量打分（cache）
	auto scoreAll = [&](bool parallel) {
		#pragma omp parallel for schedule(dynamic) num_threads(P) if(parallel)
		for (int k = 0; k < K; ++k) {
			int const v = S.internal[k];
			if (S.c0[v] < 0) continue;   // 已不是内部节点（SPR 中间态不会发生，保险）
			Scratch& sc = scratch[omp_get_thread_num()];
			int const a = S.c0[v], b = S.c1[v];
			u64 const* const bA = S.bits.data() + (size_t)a * S.W;
			u64 const* const bB = S.bits.data() + (size_t)b * S.W;
			int const sC = (int)nTaxa - S.sz[a] - S.sz[b];
			for (size_t w = 0; w < S.W; ++w) sc.t1[w] = ~(bA[w] | bB[w]);
			w3::Scratch3& sc3 = s3[omp_get_thread_num()];
			for (int g = 0; g < T; ++g)
				cache[(size_t)k * T + g] = w3::Wnode3(gv[(size_t)g], bA, bB, sc.t1.data(), S.W, sc3);
		}
	};
	scoreAll(true);
	i64 total = 0;
	for (i64 x : cache) total += x;
	total /= 2;
	std::cerr << "[spr] initial score = " << total << std::endl;

	std::mt19937_64 rng((u64)args.seed);
	vector<u64> tmpA(S.W);
	vector<i64> newW(T);

	vector<char> mark((size_t)S.n, 0);
	vector<int> allGenes((size_t)T), screenGenes;
	for (int i = 0; i < T; ++i) allGenes[(size_t)i] = i;
	if (args.screen > 0 && args.screen < T) {
		int const stride = std::max(1, T / args.screen);
		for (int i = 0; i < T; i += stride) screenGenes.push_back(i);
	}

	int kickLeft = args.kick;
	// 调试：bits 结构不变量（bits[z] == bits[c0] | bits[c1]，sz 一致）
	auto checkBits = [&](const char* tag) {
		if (!getenv("SPR_BITS")) return;
		for (int z = 0; z < S.n; ++z) {
			if (S.c0[z] < 0 || S.c1[z] < 0) continue;
			u64 const* const a = S.bits.data() + (size_t)S.c0[z] * S.W;
			u64 const* const b = S.bits.data() + (size_t)S.c1[z] * S.W;
			u64 const* const zz = S.bits.data() + (size_t)z * S.W;
			for (size_t wi = 0; wi < S.W; ++wi) {
				if (zz[wi] != (a[wi] | b[wi])) {
					std::cerr << "[BITS:" << tag << "] node=" << z << " word=" << wi
						<< " bits=" << zz[wi] << " union=" << (a[wi] | b[wi])
						<< " c0=" << S.c0[z] << " c1=" << S.c1[z] << "\n";
					std::abort();
				}
			}
		}
	};

	// 调试：全量重算并与 cache 比对
	auto verifyAll = [&](const char* tag) {
		if (!getenv("SPR_VERIFY")) return;
		for (int k = 0; k < K; ++k) {
			int const v2 = S.internal[(size_t)k];
			int const a2 = S.c0[v2], b2 = S.c1[v2];
			if (a2 < 0 || b2 < 0) continue;
			u64 const* const bA = S.bits.data() + (size_t)a2 * S.W;
			u64 const* const bB = S.bits.data() + (size_t)b2 * S.W;
			for (int gi = 0; gi < T; ++gi) {
				Scratch& sc = scratch[omp_get_thread_num()];
				for (size_t wi = 0; wi < S.W; ++wi) sc.t1[wi] = ~(bA[wi] | bB[wi]);
				i64 const wv = w3::Wnode3(gv[(size_t)gi], bA, bB, sc.t1.data(), S.W, s3[omp_get_thread_num()]);
				if (wv != cache[(size_t)k * T + gi]) {
					std::cerr << "[VERIFY:" << tag << "] node=" << v2 << " gene=" << gi
						<< " cache=" << cache[(size_t)k * T + gi] << " recomputed=" << wv
						<< " c0=" << a2 << " c1=" << b2 << " par=" << S.par[v2] << " sz=" << S.sz[v2] << "\n";
					std::abort();
				}
			}
		}
	};

	// ===== sector 相位：叶子重插入（严格按推导实现）=====
	// 摘下叶子 x（抑制父 p）：结构 + bits + sz 一起改
	auto detachLeaf = [&](int x, int& pOut, int& gOut, int& sibOut) {
		int const p = S.par[x], g = S.par[p];
		int const sib = (S.c0[p] == x) ? S.c1[p] : S.c0[p];
		if (g < 0) { S.root = sib; S.par[sib] = -1; }
		else { if (S.c0[g] == p) S.c0[g] = sib; else S.c1[g] = sib; S.par[sib] = g; }
		u64 const* const bx = S.bits.data() + (size_t)x * S.W;
		for (int z = g; z >= 0; z = S.par[z]) {
			--S.sz[z];
			u64* const bz = S.bits.data() + (size_t)z * S.W;
			for (size_t wi = 0; wi < S.W; ++wi) bz[wi] &= ~bx[wi];
		}
		pOut = p; gOut = g; sibOut = sib;
		checkBits("detach");
	};
	auto restoreDetached = [&](int x, int p, int g, int sib) {
		if (g < 0) { S.root = p; S.par[p] = -1; }
		else { if (S.c0[g] == sib) S.c0[g] = p; else S.c1[g] = p; S.par[p] = g; }
		S.c0[p] = x; S.c1[p] = sib; S.par[x] = p; S.par[sib] = p;
		// p 的 bits/sz 必须显式重建（回滚可能把 p 的 bits 留在脏状态）
		{
			u64* const bp = S.bits.data() + (size_t)p * S.W;
			u64 const* const bx = S.bits.data() + (size_t)x * S.W;
			u64 const* const bs = S.bits.data() + (size_t)sib * S.W;
			for (size_t wi = 0; wi < S.W; ++wi) bp[wi] = bx[wi] | bs[wi];
			S.sz[p] = 1 + S.sz[sib];
		}
		u64 const* const bx = S.bits.data() + (size_t)x * S.W;
		for (int z = g; z >= 0; z = S.par[z]) {
			++S.sz[z];
			u64* const bz = S.bits.data() + (size_t)z * S.W;
			for (size_t wi = 0; wi < S.W; ++wi) bz[wi] |= bx[wi];
		}
		checkBits("restore");
	};
	// 从"已摘下"状态插到 v 上方（v<0 = 根之上）；A = {w} ∪ Path(u) ∪ Path(gOld)
	auto tryInsertDetached = [&](int x, int p, int gOld, int v, bool commit) -> i64 {
		int const u = (v < 0) ? -1 : S.par[v];
		int const vv = (v < 0) ? S.root : v;
		S.c0[p] = x; S.c1[p] = vv; S.par[x] = p; S.par[vv] = p;
		if (u < 0) { S.root = p; S.par[p] = -1; }
		else { if (S.c0[u] == vv) S.c0[u] = p; else S.c1[u] = p; S.par[p] = u; }
		{
			u64* const bw = S.bits.data() + (size_t)p * S.W;
			u64 const* const bx = S.bits.data() + (size_t)x * S.W;
			u64 const* const bv = S.bits.data() + (size_t)vv * S.W;
			for (size_t wi = 0; wi < S.W; ++wi) bw[wi] = bx[wi] | bv[wi];
			S.sz[p] = 1 + S.sz[vv];
		}
		{
			u64 const* const bx = S.bits.data() + (size_t)x * S.W;
			for (int z = u; z >= 0; z = S.par[z]) {
				u64* const bz = S.bits.data() + (size_t)z * S.W;
				for (size_t wi = 0; wi < S.W; ++wi) bz[wi] |= bx[wi];
				++S.sz[z];
			}
		}
		vector<int> affected; affected.push_back(p);
		for (int z = u; z >= 0; z = S.par[z]) affected.push_back(z);
		for (int z = gOld; z >= 0; z = S.par[z]) affected.push_back(z);
		std::sort(affected.begin(), affected.end());
		affected.erase(std::unique(affected.begin(), affected.end()), affected.end());
		i64 delta = 0;
		vector<i64> newW((size_t)affected.size() * (size_t)T);
		for (size_t zi = 0; zi < affected.size(); ++zi) {
			int const z = affected[zi];
			int const a = S.c0[z], b = S.c1[z];
			if (a < 0 || b < 0) continue;
			u64 const* const bA = S.bits.data() + (size_t)a * S.W;
			u64 const* const bB = S.bits.data() + (size_t)b * S.W;
			int const sC = (int)nTaxa - S.sz[a] - S.sz[b];
			for (int gi = 0; gi < T; ++gi) {
				Scratch& sc = scratch[omp_get_thread_num()];
				for (size_t wi = 0; wi < S.W; ++wi) sc.t1[wi] = ~(bA[wi] | bB[wi]);
				newW[zi * (size_t)T + (size_t)gi] = w3::Wnode3(gv[(size_t)gi], bA, bB, sc.t1.data(), S.W,
					s3[omp_get_thread_num()]);
			}
			int const kz = idx[(size_t)z];
			for (int gi = 0; gi < T; ++gi) delta += newW[zi * (size_t)T + (size_t)gi] - cache[(size_t)kz * T + gi];
		}
		delta /= 2;
		if (commit && delta > 0) {
			for (size_t zi = 0; zi < affected.size(); ++zi) {
				int const z = affected[zi];
				int const kz = idx[(size_t)z];
				for (int gi = 0; gi < T; ++gi) cache[(size_t)kz * T + gi] = newW[zi * (size_t)T + (size_t)gi];
			}
			return delta;
		}
		checkBits("insert");
		// 回滚到"已摘下"状态
		S.c0[p] = -1; S.c1[p] = -1;
		if (u < 0) S.root = vv;
		else { if (S.c0[u] == p) S.c0[u] = vv; else S.c1[u] = vv; }
		S.par[vv] = u;
		{
			u64 const* const bx = S.bits.data() + (size_t)x * S.W;
			for (int z = u; z >= 0; z = S.par[z]) {
				u64* const bz = S.bits.data() + (size_t)z * S.W;
				for (size_t wi = 0; wi < S.W; ++wi) bz[wi] &= ~bx[wi];
				--S.sz[z];
			}
		}
		checkBits("revert");
		return delta;
	};
	// 自检 2：detach -> tryInsert(revert) -> restore 是否保持 bits 不变量
	if (getenv("SELFTEST_INSERT")) {
		int n2 = 0;
		for (int u = 0; u < S.n; ++u) {
			if (S.taxon[u] < 0) continue;
			int p, g, sib;
			detachLeaf(u, p, g, sib);
			for (int v2 = 0; v2 < S.n; ++v2) {
				if (S.taxon[v2] >= 0 && S.sz[v2] != 1) continue;
				if (v2 == u || v2 == p) continue;
				if (getenv("SELFTEST_TRACE")) std::cerr << "[trace] x=" << u << " p=" << p << " g=" << g << " v=" << v2 << "\n";
				tryInsertDetached(u, p, g, v2, false);
				checkBits("revert");
				++n2;
			}
			restoreDetached(u, p, g, sib);
			checkBits("restore2");
		}
		std::cerr << "[selftest] insert/revert " << n2 << " times: OK\n";
		return 0;
	}
	// 自检：detach/restore 是否保持 bits 不变量
	if (getenv("SELFTEST_DETACH")) {
		int bad = 0;
		for (int u = 0; u < S.n; ++u) {
			if (S.taxon[u] < 0) continue;
			int p, g, sib;
			detachLeaf(u, p, g, sib);
			checkBits("detach");
			restoreDetached(u, p, g, sib);
			checkBits("restore");
			++bad;
		}
		std::cerr << "[selftest] detach/restore on " << bad << " leaves: OK\n";
		return 0;
	}
	// 调试：--probe "xName:vName" 打印把叶子 x 移到 v 上方的精确 Δ
	if (!args.probe.empty()) {
		size_t const c = args.probe.find(':');
		string const xn = args.probe.substr(0, c), vn = args.probe.substr(c + 1);
		int xid = -1, vid = -1;
		for (int t = 0; t < (int)nTaxa; ++t) {
			string const nm = common::taxonName2ID[(size_t)t];
			if (nm == xn) xid = t;
			if (nm == vn) vid = t;
		}
		int xn2 = -1, vn2 = -1;
		for (int u = 0; u < S.n; ++u) {
			if (S.taxon[u] == xid) xn2 = u;
			if (S.taxon[u] == vid) vn2 = u;
		}
		if (xn2 < 0 || vn2 < 0) { std::cerr << "probe: taxon not found\n"; return 1; }
		int p, g, sib;
		detachLeaf(xn2, p, g, sib);
		i64 const d = tryInsertDetached(xn2, p, g, vn2, false);
		std::cout << d << "\n";
		return 0;
	}
	vector<std::pair<i64,int>> pool;
	for (int round = 1; round <= args.rounds; ++round) {
		i64 gain = 0;
		int accepted = 0;
		if (args.sectorK > 0) {
			int done = 0, acc2 = 0; i64 gain2 = 0;
			for (int kk = 0; kk < args.sectorK; ++kk) {
				int x = -1;
				for (int attempt = 0; attempt < 200; ++attempt) {
					int const cand = (int)(rng() % (u64)S.n);
					if (S.taxon[cand] >= 0) { x = cand; break; }
				}
				if (x < 0) continue;
				int p, g, sib;
				detachLeaf(x, p, g, sib);
				vector<i64> sd; i64 rootDelta = 0; vector<int> dummy;
				sector::computeDeltas(S, gv, args.sectorGenes, x, sd, rootDelta, dummy, (int)nTaxa, s3[0]);
				vector<std::pair<i64,int>> cands;
				for (int u2 = 0; u2 < S.n; ++u2) {
					if (u2 == S.root || S.taxon[u2] >= 0) continue;
					if (sd[(size_t)u2] > 0) cands.emplace_back(sd[(size_t)u2], u2);
				}
				if (rootDelta > 0) cands.emplace_back(rootDelta, -1);
				std::sort(cands.begin(), cands.end(), std::greater<>());
				if (cands.size() > 3) cands.resize(3);
				i64 best = 0; int bestV = -2;
				for (auto const& cc : cands) {
					i64 const d = tryInsertDetached(x, p, g, cc.second, false);
					if (d > best) { best = d; bestV = cc.second; }
				}
				++done;
				if (bestV != -2 && best > 0) {
					if (getenv("SPR_VERIFY")) {
						std::cerr << "[DBG] move x=" << x << " p=" << p << " g=" << g << " v=" << bestV
							<< " | Path(g):"; for (int z = g; z >= 0; z = S.par[z]) std::cerr << " " << z;
						int const uu = (bestV < 0) ? -1 : S.par[bestV];
						std::cerr << " | Path(u):"; for (int z = uu; z >= 0; z = S.par[z]) std::cerr << " " << z;
						std::cerr << std::endl;
					}
					tryInsertDetached(x, p, g, bestV, true);
					verifyAll("commit");
					total += best; gain += best; ++accepted; ++acc2; gain2 += best;
				} else {
					restoreDetached(x, p, g, sib);
				}
			}
			std::cerr << "[spr] sector: " << done << " taxa, accepted " << acc2 << " gain " << gain2 << std::endl;
		}
		if (args.target > 0.0) {
			pool.clear();
			for (int k = 0; k < K; ++k) {
				int const v = S.internal[(size_t)k];
				if (S.sz[v] < 2 || S.sz[v] > (int)nTaxa - 2) continue;
				i64 c = 0;
				i64 const* const row = cache.data() + (size_t)k * T;
				for (int g = 0; g < T; ++g) c += row[g];
				pool.emplace_back(c, v);
			}
			std::sort(pool.begin(), pool.end());
			size_t const keep = std::max<size_t>(1, (size_t)((double)pool.size() * args.target));
			pool.resize(keep);
		}
		auto ts = std::chrono::steady_clock::now();
		for (int it = 0; it < args.iters; ++it) {
			if (getenv("SPR_DEBUG")) std::cerr << "it=" << it << " ";
			// ---- 1) 选子树 x（--target>0 时从低贡献池里选；扰动阶段优先大子树）
			bool const isKick = (kickLeft > 0);
			if (isKick) --kickLeft;
			int x = -1;
			if (isKick) {
				int const minSz = std::max(5, (int)nTaxa / 200);
				for (int attempt = 0; attempt < 200; ++attempt) {
					int const cand = S.internal[(size_t)(rng() % (u64)K)];
					if (S.sz[cand] >= minSz && S.sz[cand] <= (int)nTaxa - minSz) { x = cand; break; }
				}
			}
			if (x < 0) {
				if (args.target > 0.0 && !pool.empty()) x = pool[(size_t)(rng() % (u64)pool.size())].second;
				else if (args.allowLeaf) x = (int)(rng() % (u64)S.n);
				else x = S.internal[(size_t)(rng() % (u64)K)];
			}
			if (S.sz[x] < (args.allowLeaf ? 1 : 2) || S.sz[x] > (int)nTaxa - 2) continue;
			int const p = S.par[x];
			if (p < 0) continue;
			int const g = S.par[p];
			int const sib = (S.c0[p] == x) ? S.c1[p] : S.c0[p];

			// ---- 2) 保存受影响节点的 bits/sz（旧路径 + 新路径，去重；恢复时按保存值还原）
			vector<int> path1;
			for (int z = g; z >= 0; z = S.par[z]) path1.push_back(z);
			vector<int> savedNodes;
			vector<u64> savedBits;
			vector<int> savedSz;
			vector<char> saveMark((size_t)S.n, 0);
			vector<int> saveIdx((size_t)S.n, -1);
			auto saveNode = [&](int z) {
				if (z < 0 || saveMark[(size_t)z]) return;
				saveMark[(size_t)z] = 1;
				saveIdx[(size_t)z] = (int)savedNodes.size();
				size_t const at = savedBits.size();
				savedBits.resize(at + S.W);
				memcpy(savedBits.data() + at, S.bits.data() + (size_t)z * S.W, S.W * sizeof(u64));
				savedNodes.push_back(z);
				savedSz.push_back(S.sz[z]);
			};
			for (int z : path1) saveNode(z);
			saveNode(p);   // p 会被复用为 w，回滚时要恢复它的 bits/sz
			int const c0p = S.c0[p], c1p = S.c1[p], parP = S.par[p];
			int const c0g = (g >= 0 ? S.c0[g] : -1), c1g = (g >= 0 ? S.c1[g] : -1);
			int const parSib = S.par[sib];
			int const rootOld = S.root;

			if (getenv("SPR_D8") && it == 8) {
				std::cerr << "  c0[p]=" << S.c0[p] << " c1[p]=" << S.c1[p] << " sib=" << sib
					<< " par[sib]=" << S.par[sib] << " root=" << S.root << " n=" << S.n << "\n";
				std::cerr << "  chain from 349:"; for (int q = 349, c = 0; q >= 0 && c < 12; q = S.par[q], ++c) std::cerr << " " << q;
				std::cerr << "\n  chain from 1:"; for (int q = 1, c = 0; q >= 0 && c < 12; q = S.par[q], ++c) std::cerr << " " << q;
				std::cerr << std::endl;
			}
			if (getenv("SPR_D0") && it == 0) {
				std::cerr << "  p=" << p << " c0[p]=" << S.c0[p] << " c1[p]=" << S.c1[p]
					<< " g=" << g << " c0[g]=" << (g>=0?S.c0[g]:-1) << " c1[g]=" << (g>=0?S.c1[g]:-1)
					<< " sib=" << sib << "\n";
			}
			// ---- 3) 摘下 x
			if (g < 0) { S.root = sib; S.par[sib] = -1; }
			else { if (S.c0[g] == p) S.c0[g] = sib; else S.c1[g] = sib; S.par[sib] = g; }
			{
				u64 const* const bx = S.bits.data() + (size_t)x * S.W;
				for (int z : path1) {
					u64* const bz = S.bits.data() + (size_t)z * S.W;
					for (size_t w = 0; w < S.W; ++w) bz[w] &= ~bx[w];
					S.sz[z] -= S.sz[x];
				}
			}

			// ---- 4) 选重挂点 v（不在 subtree(x) 里，且不是 p/ x）
			int v = -1;
			if (args.localDepth > 0 && g >= 0) {
				// BFS 收集 g 附近 localDepth 条边内的节点
				vector<int> nb; nb.push_back(g);
				for (int d = 0, head = 0; d < args.localDepth && head < (int)nb.size(); ++d) {
					int const tail = (int)nb.size();
					for (; head < tail; ++head) {
						int const z = nb[(size_t)head];
						if (S.c0[z] >= 0) { nb.push_back(S.c0[z]); nb.push_back(S.c1[z]); }
						if (S.par[z] >= 0) nb.push_back(S.par[z]);
					}
				}
				std::shuffle(nb.begin(), nb.end(), rng);
				for (int cand : nb) {
					if (cand == p || cand == x) continue;
					u64 const* const bc = S.bits.data() + (size_t)cand * S.W;
					u64 const* const bx = S.bits.data() + (size_t)x * S.W;
					bool inX = true;
					for (size_t w = 0; w < S.W; ++w) if ((bc[w] & bx[w]) != bc[w]) { inX = false; break; }
					if (inX) continue;
					v = cand; break;
				}
			}
			for (int attempt = 0; attempt < 80 && v < 0; ++attempt) {
				int const cand = (int)(rng() % (u64)S.n);
				if (cand == p || cand == x) continue;
				u64 const* const bc = S.bits.data() + (size_t)cand * S.W;
				u64 const* const bx = S.bits.data() + (size_t)x * S.W;
				bool inX = true;
				for (size_t w = 0; w < S.W; ++w) if ((bc[w] & bx[w]) != bc[w]) { inX = false; break; }
				if (inX) continue;
				v = cand; break;
			}
			if (v < 0) {   // 回滚
				if (g < 0) { S.root = rootOld; S.par[sib] = parSib; }
				else { S.c0[g] = c0g; S.c1[g] = c1g; S.par[sib] = parSib; }
				for (size_t i = 0; i < savedNodes.size(); ++i) {
					int const z = savedNodes[i];
					memcpy(S.bits.data() + (size_t)z * S.W, savedBits.data() + i * S.W, S.W * sizeof(u64));
					S.sz[z] = savedSz[i];
				}
				continue;
			}
			int const u = S.par[v];
			if (getenv("SPR_D55") && it == 55) {
				std::cerr << "  BEFORE: p=" << p << " c0[p]=" << S.c0[p] << " c1[p]=" << S.c1[p]
					<< " g=" << g << " c0[g]=" << (g>=0?S.c0[g]:-1) << " c1[g]=" << (g>=0?S.c1[g]:-1)
					<< " sib=" << sib << " v=" << v << " u=" << u
					<< " c0[u]=" << (u>=0?S.c0[u]:-1) << " c1[u]=" << (u>=0?S.c1[u]:-1) << "\n";
			}
			if (getenv("SPR_DEBUG")) std::cerr << "x=" << x << " p=" << p << " g=" << g << " v=" << v << " u=" << u << "\n";
			int const c0u = (u >= 0 ? S.c0[u] : -1), c1u = (u >= 0 ? S.c1[u] : -1), parV = S.par[v];

			// ---- 5) 重挂：w 复用 p 的槽位
			int const w0 = p;
			S.c0[w0] = x; S.c1[w0] = v; S.par[x] = w0; S.par[v] = w0;
			if (u < 0) { S.root = w0; S.par[w0] = -1; }
			else { if (S.c0[u] == v) S.c0[u] = w0; else S.c1[u] = w0; S.par[w0] = u; }
			{
				u64* const bw = S.bits.data() + (size_t)w0 * S.W;
				u64 const* const bx = S.bits.data() + (size_t)x * S.W;
				u64 const* const bv = S.bits.data() + (size_t)v * S.W;
				for (size_t wi = 0; wi < S.W; ++wi) bw[wi] = bx[wi] | bv[wi];
				S.sz[w0] = S.sz[x] + S.sz[v];
			}
			vector<int> path2;
			for (int z = u; z >= 0; z = S.par[z]) path2.push_back(z);
			for (int z : path2) saveNode(z);   // 保存"加 x 之前"的状态（path1 节点已存过原值）
			{
				u64 const* const bx = S.bits.data() + (size_t)x * S.W;
				for (int z : path2) {
					u64* const bz = S.bits.data() + (size_t)z * S.W;
					for (size_t wi = 0; wi < S.W; ++wi) bz[wi] |= bx[wi];
					S.sz[z] += S.sz[x];
				}
			}

			// ---- 6) 受影响内部节点 = 两条路径 ∪ {w}，去重
			vector<int> affected;
			auto push = [&](int z) {
				if (z < 0 || idx[(size_t)z] < 0) return;
				if (mark[(size_t)z]) return;
				mark[(size_t)z] = 1;
				affected.push_back(z);
			};
			for (int z : path1) push(z);
			for (int z : path2) push(z);
			push(w0);

			for (int z : affected) mark[(size_t)z] = 0;   // 去重完成，清标记
			// ---- 7) 评估：单个并行区处理 (受影响节点 × 基因)
			std::vector<i64> fullW;
			auto evalMove = [&](vector<int> const& glist, std::vector<i64>* outW) -> i64 {
				size_t const nn = affected.size(), GG = glist.size();
				if (outW) outW->assign(nn * GG, 0);
				i64 d = 0;
				long long const total_pairs = (long long)(nn * GG);
				#pragma omp parallel for schedule(dynamic) num_threads(P) reduction(+:d)
				for (long long t = 0; t < total_pairs; ++t) {
					size_t const zi = (size_t)t / GG, gi = (size_t)t % GG;
					int const z = affected[zi];
					int const a = S.c0[z], b = S.c1[z];
					if (a < 0 || b < 0) continue;
					u64 const* const bA = S.bits.data() + (size_t)a * S.W;
					u64 const* const bB = S.bits.data() + (size_t)b * S.W;
					int const sC = (int)nTaxa - S.sz[a] - S.sz[b];
					Scratch& sc = scratch[omp_get_thread_num()];
					for (size_t wi = 0; wi < S.W; ++wi) sc.t1[wi] = ~(bA[wi] | bB[wi]);
					i64 const wv = w3::Wnode3(gv[(size_t)glist[gi]], bA, bB, sc.t1.data(), S.W,
						s3[omp_get_thread_num()]);
					if (outW) (*outW)[(size_t)t] = wv;
					d += wv - cache[(size_t)idx[(size_t)z] * T + glist[gi]];
				}
				return d;
			};
			i64 delta = 0;
			bool pass = true;
			if (!screenGenes.empty()) {
				i64 const ds = evalMove(screenGenes, nullptr);
				if (ds <= 0) pass = false;
			}
			if (pass) delta = evalMove(allGenes, &fullW);
			delta /= 2;

			if (delta > 0 || isKick) {
				if (!fullW.empty()) {
					size_t const GG = allGenes.size();
					for (size_t zi = 0; zi < affected.size(); ++zi) {
						int const kz = idx[(size_t)affected[zi]];
						if (kz < 0) continue;
						for (size_t gi = 0; gi < GG; ++gi)
							cache[(size_t)kz * T + allGenes[gi]] = fullW[zi * GG + gi];
					}
				}
				total += delta;
				gain += delta;
				++accepted;
			}
			else {
				if (g < 0) { S.root = rootOld; S.par[sib] = parSib; }
				else { S.c0[g] = c0g; S.c1[g] = c1g; S.par[sib] = parSib; }
				for (size_t i = 0; i < savedNodes.size(); ++i) {
					int const z = savedNodes[i];
					memcpy(S.bits.data() + (size_t)z * S.W, savedBits.data() + i * S.W, S.W * sizeof(u64));
					S.sz[z] = savedSz[i];
				}
				S.c0[p] = c0p; S.c1[p] = c1p; S.par[p] = parP;
				S.par[x] = p;
				if (v != sib) S.par[v] = parV;   // v==sib 时 par[sib] 已按原值恢复
				if (u < 0) S.root = rootOld;
				else if (u != g) { S.c0[u] = c0u; S.c1[u] = c1u; }   // u==g 时上面已按原值恢复
			}
			for (int z : affected) mark[(size_t)z] = 0;
			if (getenv("SPR_BITS")) {
				for (int z = 0; z < S.n; ++z) {
					if (S.c0[z] < 0 || S.c1[z] < 0) continue;
					u64* zb = S.bits.data() + (size_t)z * S.W;
					u64 const* a = S.bits.data() + (size_t)S.c0[z] * S.W;
					u64 const* b = S.bits.data() + (size_t)S.c1[z] * S.W;
					int sum = 0;
					for (size_t w2 = 0; w2 < S.W; ++w2) {
						if (zb[w2] != (a[w2] | b[w2])) {
							std::cerr << "[spr] BITS MISMATCH at it=" << it << " node=" << z
								<< " c0=" << S.c0[z] << " c1=" << S.c1[z] << " par=" << S.par[z]
								<< " sz=" << S.sz[z] << "\n";
							std::cerr << "   move: x=" << x << " p=" << p << " g=" << g << " v=" << v << " u=" << u << " w0=" << w0 << "\n";
							std::abort();
						}
						sum += __builtin_popcountll(a[w2]) + __builtin_popcountll(b[w2]);
					}
					if (sum != S.sz[z]) {
						std::cerr << "[spr] SZ MISMATCH at it=" << it << " node=" << z << " sz=" << S.sz[z] << " sum=" << sum << "\n";
						std::abort();
					}
				}
			}
			if (getenv("SPR_VALID") && !validateTree(S)) {
				std::cerr << "[spr] INVALID TREE after it=" << it << " x=" << x << " p=" << p << " g=" << g
					<< " v=" << v << " u=" << u << " delta=" << delta << std::endl;
				return 3;
			}
		}
		std::cerr << "[spr] round " << round << ": accepted=" << accepted << " gain=" << gain
			<< " total=" << total << "  ("
			<< std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - ts).count()
			<< " s)" << std::endl;
		if (!args.output.empty()) {
			std::ofstream fout(args.output);
			fout << toNewick(S, S.root) << ";\n";
		}
		if (accepted == 0) break;
	}

	string const nw = toNewick(S, S.root) + ";";
	if (args.output.empty()) std::cout << nw << std::endl;
	else { std::ofstream fout(args.output); fout << nw << "\n"; }
	std::cerr << "[spr] final score = " << total << std::endl;
	return 0;
}
