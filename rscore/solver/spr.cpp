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

struct Args {
	string species, input, output;
	int threads = 8, rounds = 50, iters = 2000;
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
			for (int g = 0; g < T; ++g)
				cache[(size_t)k * T + g] = Wnode(genes[g], bA, bB, sc.t1.data(),
					S.sz[a], S.sz[b], sC, S.W, sc);
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

	for (int round = 1; round <= args.rounds; ++round) {
		i64 gain = 0;
		int accepted = 0;
		auto ts = std::chrono::steady_clock::now();
		for (int it = 0; it < args.iters; ++it) {
			if (getenv("SPR_DEBUG")) std::cerr << "it=" << it << " ";
			// ---- 1) 选子树 x
			int const kx = (int)(rng() % (u64)K);
			int const x = S.internal[(size_t)kx];
			if (S.sz[x] < 2 || S.sz[x] > (int)nTaxa - 2) continue;
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
			auto saveNode = [&](int z) {
				if (z < 0 || saveMark[(size_t)z]) return;
				saveMark[(size_t)z] = 1;
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
			for (int attempt = 0; attempt < 80; ++attempt) {
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

			// ---- 7) 重算 W 并累加 delta
			i64 delta = 0;
			for (int z : affected) {
				int const kz = idx[(size_t)z];
				int const a = S.c0[z], b = S.c1[z];
				if (a < 0 || b < 0) continue;
				u64 const* const bA = S.bits.data() + (size_t)a * S.W;
				u64 const* const bB = S.bits.data() + (size_t)b * S.W;
				int const sC = (int)nTaxa - S.sz[a] - S.sz[b];
				#pragma omp parallel for schedule(dynamic) num_threads(P)
				for (int gi = 0; gi < T; ++gi) {
					Scratch& sc = scratch[omp_get_thread_num()];
					for (size_t wi = 0; wi < S.W; ++wi) sc.t1[wi] = ~(bA[wi] | bB[wi]);
					newW[gi] = Wnode(genes[gi], bA, bB, sc.t1.data(), S.sz[a], S.sz[b], sC, S.W, sc);
				}
				for (int gi = 0; gi < T; ++gi) delta += newW[gi] - cache[(size_t)kz * T + gi];
			}
			delta /= 2;

			if (delta > 0) {
				for (int z : affected) {
					int const kz = idx[(size_t)z];
					int const a = S.c0[z], b = S.c1[z];
					if (a < 0 || b < 0) continue;
					u64 const* const bA = S.bits.data() + (size_t)a * S.W;
					u64 const* const bB = S.bits.data() + (size_t)b * S.W;
					int const sC = (int)nTaxa - S.sz[a] - S.sz[b];
					#pragma omp parallel for schedule(dynamic) num_threads(P)
					for (int gi = 0; gi < T; ++gi) {
						Scratch& sc = scratch[omp_get_thread_num()];
						for (size_t wi = 0; wi < S.W; ++wi) sc.t1[wi] = ~(bA[wi] | bB[wi]);
						cache[(size_t)kz * T + gi] = Wnode(genes[gi], bA, bB, sc.t1.data(),
							S.sz[a], S.sz[b], sC, S.W, sc);
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
