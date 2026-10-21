// recombine: 把多棵局部最优物种树的三分划并集做成 DAG，用 DP 找最优兼容组合。
//
// 目标：score = Σ_v W(A_v,B_v,C_v)/2。每棵树贡献若干 (cluster) 及其 split；
// 候选 split 只取"在某棵输入树里真实出现过"的，保证 DAG 规模可控。
//
// 用法: ./recombine -i genes.trees -s tree1,tree2,... -o out.tree [-t N] [--genes K]

#include "newick.hpp"
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <fstream>
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

struct InTree {
	int n = 0, root = -1;
	vector<int> par, c0, c1, taxon;
};

} // namespace

int main(int argc, char** argv) {
	string genesFile, speciesFiles, outFile;
	int threads = 8, Tgenes = 100;
	for (int i = 1; i < argc; ++i) {
		string s = argv[i];
		auto nx = [&]() -> string { if (i + 1 < argc) return argv[++i]; std::exit(2); };
		if (s == "-i") genesFile = nx();
		else if (s == "-s") speciesFiles = nx();
		else if (s == "-o") outFile = nx();
		else if (s == "-t") threads = std::stoi(nx());
		else if (s == "--genes") Tgenes = std::stoi(nx());
	}
	if (genesFile.empty() || speciesFiles.empty()) {
		std::cerr << "usage: recombine -i genes.trees -s t1,t2,... -o out.tree [-t N] [--genes K]" << std::endl;
		return 2;
	}

	// 基因树
	vector<GeneTree> genes;
	{
		std::ifstream fin(genesFile);
		if (!fin.is_open()) { std::cerr << "Error: cannot open " << genesFile << std::endl; return 2; }
		string line;
		while (std::getline(fin, line)) {
			if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
			newick::Newick nw(line);
			vector<newick::Newick::Node const*> order, st{ nw.root() };
			while (!st.empty()) { auto* u = st.back(); st.pop_back(); order.push_back(u);
				for (auto it = u->children.rbegin(); it != u->children.rend(); ++it) st.push_back(*it); }
			std::unordered_map<newick::Newick::Node const*, int> id;
			for (size_t i = 0; i < order.size(); ++i) id[order[i]] = (int)i;
			GeneTree g; g.n = (int)order.size();
			g.start.assign(g.n + 1, 0); g.solid.assign(g.n, 1); g.taxon.assign(g.n, -1);
			for (size_t i = 0; i < order.size(); ++i) g.start[i+1] = g.start[i] + (int)order[i]->children.size();
			for (size_t i = 0; i < order.size(); ++i) {
				auto* u = order[i];
				if (u->children.empty()) g.taxon[i] = (int)common::taxonName2ID[u->name];
				else for (auto* c : u->children) { int cid = id[c]; g.kids.push_back(cid); g.solid[cid] = (c->name == "D") ? 0 : 1; }
			}
			g.root = 0;
			genes.push_back(g);
		}
	}
	if (genes.empty()) { std::cerr << "Error: no gene trees" << std::endl; return 2; }
	size_t const nTaxa = common::taxonName2ID.nTaxa();
	size_t const W = (nTaxa + 63) / 64;
	int const T = std::min<int>(Tgenes, (int)genes.size());

	// 物种树
	vector<InTree> trees;
	{
		std::ifstream fin(speciesFiles);
		if (!fin.is_open()) { std::cerr << "Error: cannot open " << speciesFiles << std::endl; return 2; }
		string line;
		while (std::getline(fin, line)) {
			if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
			newick::Newick nw(line);
			vector<newick::Newick::Node const*> order, st{ nw.root() };
			while (!st.empty()) { auto* u = st.back(); st.pop_back(); order.push_back(u);
				for (auto it = u->children.rbegin(); it != u->children.rend(); ++it) st.push_back(*it); }
			std::unordered_map<newick::Newick::Node const*, int> id;
			for (size_t i = 0; i < order.size(); ++i) id[order[i]] = (int)i;
			InTree t; t.n = (int)order.size(); t.par.assign(t.n, -1); t.c0.assign(t.n, -1); t.c1.assign(t.n, -1); t.taxon.assign(t.n, -1);
			for (size_t i = 0; i < order.size(); ++i) {
				auto* u = order[i];
				if (u->children.empty()) t.taxon[i] = (int)common::taxonName2ID[u->name];
				else for (auto* c : u->children) {
					int cid = id[c]; t.par[cid] = (int)i;
					if (t.c0[i] < 0) t.c0[i] = cid; else t.c1[i] = cid;
				}
			}
			t.root = 0;
			trees.push_back(t);
		}
	}
	if (trees.empty()) { std::cerr << "Error: no species trees" << std::endl; return 2; }

	// 收集 cluster（bitset）与 split
	vector<vector<u64>> clusters;
	std::unordered_map<string, int> clusterId;
	vector<int> clusterSize;
	vector<vector<std::pair<int,int>>> splits;   // cluster -> [(child1, child2)]
	vector<i64> clusterWsum;                       // 该 cluster 作为节点时 ΣW（填进 clusterWsum）
	auto keyOf = [&](u64 const* b) { return string((char const*)b, W * sizeof(u64)); };
	auto getCluster = [&](u64 const* b, int sz) -> int {
		string k = keyOf(b);
		auto it = clusterId.find(k);
		if (it != clusterId.end()) return it->second;
		int id = (int)clusters.size();
		clusters.emplace_back(b, b + W);
		clusterSize.push_back(sz);
		splits.emplace_back();
		clusterWsum.push_back(0);
		clusterId.emplace(std::move(k), id);
		return id;
	};

	vector<vector<u64>> treeBits(trees.size());
	for (size_t ti = 0; ti < trees.size(); ++ti) {
		InTree const& t = trees[ti];
		treeBits[ti].assign((size_t)t.n * W, 0);
		vector<int> sz(t.n, 0);
		for (int u = t.n - 1; u >= 0; --u) {
			if (t.taxon[u] >= 0) { sz[u] = 1; treeBits[ti][(size_t)u*W + ((size_t)t.taxon[u] >> 6)] |= 1ULL << (t.taxon[u] & 63); }
			else {
				for (size_t k = 0; k < W; ++k) treeBits[ti][(size_t)u*W + k] = treeBits[ti][(size_t)t.c0[u]*W + k] | treeBits[ti][(size_t)t.c1[u]*W + k];
				sz[u] = sz[t.c0[u]] + sz[t.c1[u]];
			}
		}
		for (int u = 0; u < t.n; ++u) {
			if (t.c0[u] < 0) continue;
			u64 const* b = treeBits[ti].data() + (size_t)u * W;
			int cid = getCluster(b, sz[u]);
			int a = getCluster(treeBits[ti].data() + (size_t)t.c0[u]*W, sz[t.c0[u]]);
			int c = getCluster(treeBits[ti].data() + (size_t)t.c1[u]*W, sz[t.c1[u]]);
			auto& sp = splits[cid];
			bool seen = false;
			for (auto& p : sp) if ((p.first == a && p.second == c) || (p.first == c && p.second == a)) { seen = true; break; }
			if (!seen) sp.emplace_back(a, c);
		}
	}
	std::cerr << "[recombine] clusters=" << clusters.size() << " taxa=" << nTaxa << std::endl;

	// 按大小排序做 DP
	vector<int> order(clusters.size());
	for (size_t i = 0; i < order.size(); ++i) order[i] = (int)i;
	std::sort(order.begin(), order.end(), [&](int x, int y) { return clusterSize[x] < clusterSize[y]; });

	int const fullId = [&]() {
		vector<u64> b(W, 0);
		for (size_t i = 0; i < nTaxa; ++i) b[i >> 6] |= 1ULL << (i & 63);
		return getCluster(b.data(), (int)nTaxa);
	}();

	vector<Scratch> scratch(threads);
	int maxN = 1; for (auto const& g : genes) maxN = std::max(maxN, g.n);
	for (auto& s : scratch) s.resize(maxN, W);

	vector<i64> best(clusters.size(), LLONG_MIN);
	for (size_t i = 0; i < clusters.size(); ++i)
		if (clusterSize[i] == 1) best[i] = 0;   // 叶子：基态
	vector<int> bestSplit(clusters.size(), -1);
	{
		Scratch& sc = scratch[0];
		vector<u64> cBits(W);
		for (int id : order) {
			auto const& sp = splits[id];
			if (id == fullId) {
				// root：C 为空
				for (size_t si = 0; si < sp.size(); ++si) {
					int const a = sp[si].first, b = sp[si].second;
					if (best[a] == LLONG_MIN || best[b] == LLONG_MIN) continue;
					i64 w = 0;
					for (size_t k = 0; k < W; ++k) cBits[k] = 0;
					for (int g = 0; g < T; ++g)
						w += Wnode(genes[g], clusters[a].data(), clusters[b].data(), cBits.data(),
							clusterSize[a], clusterSize[b], 0, W, sc);
					i64 val = best[a] + best[b] + w;
					if (val > best[id]) { best[id] = val; bestSplit[id] = (int)si; }
				}
				continue;
			}
			for (size_t si = 0; si < sp.size(); ++si) {
				int const a = sp[si].first, b = sp[si].second;
				if (best[a] == LLONG_MIN || best[b] == LLONG_MIN) continue;
				int const sC = (int)nTaxa - clusterSize[a] - clusterSize[b];
				for (size_t k = 0; k < W; ++k) cBits[k] = ~(clusters[a][k] | clusters[b][k]);
				i64 w = 0;
				for (int g = 0; g < T; ++g)
					w += Wnode(genes[g], clusters[a].data(), clusters[b].data(), cBits.data(),
						clusterSize[a], clusterSize[b], sC, W, sc);
				i64 val = best[a] + best[b] + w;
				if (val > best[id]) { best[id] = val; bestSplit[id] = (int)si; }
			}
		}
	}
	std::cerr << "[recombine] DP score = " << best[fullId] / 2 << std::endl;

	// 输出
	std::function<string(int)> emit = [&](int id) -> string {
		if (splits[id].empty()) {
			return "";   // 叶子：由调用者处理
		}
		int si = bestSplit[id];
		int a = splits[id][si].first, b = splits[id][si].second;
		auto sub = [&](int x) -> string {
			if (splits[x].empty() && clusterSize[x] == 1) {
				for (size_t i = 0; i < nTaxa; ++i) if ((clusters[x][i >> 6] >> (i & 63)) & 1ULL) return common::taxonName2ID[i];
			}
			return emit(x);
		};
		return "(" + sub(a) + "," + sub(b) + ")";
	};
	string nw = emit(fullId) + ";";
	if (outFile.empty()) std::cout << nw << std::endl;
	else { std::ofstream f(outFile); f << nw << "\n"; }
	return 0;
}
