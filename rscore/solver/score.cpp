// score: 仅评分模式（助教用）。读入一棵估计物种树 + 一个基因树文件，输出可靠四叶树分数（JSON）。
//
// 用法:
//   ./score -t species.tree -g genes.trees[.gz] [-j THREADS] [--quiet]
//
// 输出（stdout，JSON）:
//   {"score":N,"n_taxa":..,"n_genes":..,"n_scored_genes":..,"n_skipped_genes":..,
//    "n_unused_taxa":..,"n_ignored_taxa":..,"runtime_sec":..,"per_gene":[...]}
//
// 退出码: 0 成功；1 输入非法（物种树非二叉/重名/空文件等）；2 参数错误。

#include "newick.hpp"
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
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
	vector<i64> c1, a1, m1, G1, H1, c2, a2, m2, G2, H2;
	vector<u64> t1, t2;
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
	if (sA >= 2 && sB + sC >= 2) { for (size_t w = 0; w < W; ++w) t[w] = B[w] | C[w]; gA = (i64)evalF(g, A, t, sc); }
	if (sB >= 2 && sA + sC >= 2) { for (size_t w = 0; w < W; ++w) t[w] = A[w] | C[w]; gB = (i64)evalF(g, B, t, sc); }
	if (sC >= 2 && sA + sB >= 2) { for (size_t w = 0; w < W; ++w) t[w] = A[w] | B[w]; gC = (i64)evalF(g, C, t, sc); }
	return gA + gB + gC - 2 * (fAB + fAC + fBC);
}

struct SpeciesTree {
	int n = 0, root = 0;
	vector<int> par, c0, c1, taxon, sz;
	vector<int> internal;
	vector<u64> bits;
	size_t W = 0;
};

// 返回空字符串表示成功，否则为错误信息。
string buildSpecies(newick::Newick::Node const* const root, size_t nTaxa, SpeciesTree& t) {
	vector<newick::Newick::Node const*> order, st{ root };
	while (!st.empty()) { auto* u = st.back(); st.pop_back(); order.push_back(u);
		for (auto it = u->children.rbegin(); it != u->children.rend(); ++it) st.push_back(*it); }
	std::unordered_map<newick::Newick::Node const*, int> id;
	id.reserve(order.size() * 2);
	for (size_t i = 0; i < order.size(); ++i) id[order[i]] = (int)i;
	t.n = (int)order.size();
	t.par.assign(t.n, -1); t.c0.assign(t.n, -1); t.c1.assign(t.n, -1); t.taxon.assign(t.n, -1);
	t.root = 0;
	t.W = (nTaxa + 63) / 64;
	t.bits.assign((size_t)t.n * t.W, 0);
	t.sz.assign(t.n, 0);
	std::set<string> seen;
	for (size_t i = 0; i < order.size(); ++i) {
		auto* u = order[i];
		if (u->children.empty()) {
			if (u->name.empty()) return "species tree has an empty leaf label";
			if (!seen.insert(u->name).second) return "duplicate leaf label in species tree: " + u->name;
			t.taxon[i] = (int)common::taxonName2ID[u->name];
		} else {
			if (u->children.size() != 2)
				return "species tree must be binary (node with " + std::to_string(u->children.size())
					+ " children); polytomous species trees are not supported";
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
		} else {
			u64 const* const a = t.bits.data() + (size_t)t.c0[i] * t.W;
			u64 const* const b = t.bits.data() + (size_t)t.c1[i] * t.W;
			for (size_t w = 0; w < t.W; ++w) dst[w] = a[w] | b[w];
			t.sz[i] = t.sz[t.c0[i]] + t.sz[t.c1[i]];
		}
	}
	return {};
}

} // namespace

int main(int argc, char** argv) {
	string treeFile, genesFile;
	int threads = 0;
	bool quiet = false;
	for (int i = 1; i < argc; ++i) {
		string const s = argv[i];
		auto nx = [&]() -> string { if (i + 1 < argc) return argv[++i]; std::exit(2); };
		if (s == "-t" || s == "--tree") treeFile = nx();
		else if (s == "-g" || s == "--genes") genesFile = nx();
		else if (s == "-j" || s == "--threads") threads = std::stoi(nx());
		else if (s == "-q" || s == "--quiet") quiet = true;
		else if (s == "-h" || s == "--help") {
			std::cout << "usage: score -t species.tree -g genes.trees[.gz] [-j THREADS] [-q]\n";
			return 0;
		}
		else { std::cerr << "Error: unknown argument " << s << "\n"; return 2; }
	}
	if (treeFile.empty() || genesFile.empty()) {
		std::cerr << "Error: -t <species tree> and -g <gene trees> are required\n";
		return 2;
	}

	// 物种树：必须恰好一棵
	vector<newick::Newick> spTrees;
	{
		std::ifstream fin(treeFile);
		if (!fin.is_open()) { std::cerr << "Error: cannot open species tree file: " << treeFile << "\n"; return 1; }
		string line;
		while (std::getline(fin, line)) {
			if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
			try { spTrees.emplace_back(line); }
			catch (std::exception const& e) { std::cerr << "Error: bad Newick in " << treeFile << ": " << e.what() << "\n"; return 1; }
		}
	}
	if (spTrees.empty()) { std::cerr << "Error: no species tree found in " << treeFile << "\n"; return 1; }
	if (spTrees.size() != 1) {
		std::cerr << "Error: expected exactly one species tree in " << treeFile
			<< ", found " << spTrees.size() << "\n";
		return 1;
	}

	// 基因树
	vector<GeneTree> allGenes;
	vector<int> skipped;
	{
		std::ifstream fin(genesFile);
		if (!fin.is_open()) { std::cerr << "Error: cannot open gene tree file: " << genesFile << "\n"; return 1; }
		string line; int idx = 0;
		while (std::getline(fin, line)) {
			if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
			++idx;
			newick::Newick nw(line);
			GeneTree g = buildGene(nw.root());
			int nLeaves = 0;
			for (int u = 0; u < g.n; ++u) if (g.taxon[u] >= 0) ++nLeaves;
			if (nLeaves < 4) { skipped.push_back(idx); continue; }
			allGenes.push_back(std::move(g));
		}
	}
	if (allGenes.empty()) { std::cerr << "Error: no usable gene tree (>=4 taxa) in " << genesFile << "\n"; return 1; }

	// 先把物种树的叶子名注册进 taxon 表，nTaxa 才会包含"只出现在物种树里"的类群
	{
		vector<newick::Newick::Node const*> st{ spTrees[0].root() };
		while (!st.empty()) {
			auto* u = st.back(); st.pop_back();
			if (u->children.empty()) { if (!u->name.empty()) (void)common::taxonName2ID[u->name]; }
			else for (auto* c : u->children) st.push_back(c);
		}
	}
	size_t const nTaxa = common::taxonName2ID.nTaxa();
	SpeciesTree S;
	{
		string const err = buildSpecies(spTrees[0].root(), nTaxa, S);
		if (!err.empty()) { std::cerr << "Error: " << err << "\n"; return 1; }
	}

	// 类群覆盖检查
	vector<char> inGene(nTaxa, 0);
	for (auto const& g : allGenes)
		for (int u = 0; u < g.n; ++u)
			if (g.taxon[u] >= 0) inGene[(size_t)g.taxon[u]] = 1;
	vector<char> inSpecies(nTaxa, 0);
	for (int u = 0; u < S.n; ++u)
		if (S.taxon[u] >= 0) inSpecies[(size_t)S.taxon[u]] = 1;
	int nUnused = 0, nIgnored = 0;
	for (size_t i = 0; i < nTaxa; ++i) {
		if (inSpecies[i] && !inGene[i]) ++nUnused;
		if (!inSpecies[i] && inGene[i]) ++nIgnored;
	}
	if (!quiet) {
		if (nUnused) std::cerr << "Warning: " << nUnused << " species-tree taxa never occur in any gene tree (they contribute nothing)\n";
		if (nIgnored) std::cerr << "Warning: " << nIgnored << " gene-tree taxa are absent from the species tree (ignored)\n";
		if (!skipped.empty()) std::cerr << "Warning: skipped " << skipped.size() << " gene tree(s) with <4 taxa\n";
	}

	int const P = threads > 0 ? threads : omp_get_max_threads();
	int const T = (int)allGenes.size();
	vector<Scratch> scratch((size_t)P);
	int maxN = 1; for (auto const& g : allGenes) maxN = std::max(maxN, g.n);
	for (auto& sc : scratch) sc.resize(maxN, S.W);

	auto const t0 = std::chrono::steady_clock::now();
	vector<i64> perGene((size_t)T, 0);
	#pragma omp parallel for schedule(dynamic) num_threads(P)
	for (int gi = 0; gi < T; ++gi) {
		Scratch& sc = scratch[(size_t)omp_get_thread_num()];
		GeneTree const& g = allGenes[(size_t)gi];
		i64 sum = 0;
		for (int v : S.internal) {
			int const a = S.c0[v], b = S.c1[v];
			u64 const* const bA = S.bits.data() + (size_t)a * S.W;
			u64 const* const bB = S.bits.data() + (size_t)b * S.W;
			int const sC = (int)nTaxa - S.sz[a] - S.sz[b];
			for (size_t w = 0; w < S.W; ++w) sc.t1[w] = ~(bA[w] | bB[w]);
			sum += Wnode(g, bA, bB, sc.t1.data(), S.sz[a], S.sz[b], sC, S.W, sc);
		}
		perGene[(size_t)gi] = sum / 2;
	}
	i64 score = 0;
	for (i64 x : perGene) score += x;
	double const secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

	// JSON 输出
	std::cout << "{\"score\":" << score
		<< ",\"n_taxa\":" << (int)nTaxa
		<< ",\"n_genes\":" << T
		<< ",\"n_scored_genes\":" << T
		<< ",\"n_skipped_genes\":" << (int)skipped.size()
		<< ",\"n_unused_taxa\":" << nUnused
		<< ",\"n_ignored_taxa\":" << nIgnored
		<< ",\"runtime_sec\":" << std::fixed << std::setprecision(3) << secs
		<< ",\"per_gene\":[";
	for (int i = 0; i < T; ++i) { if (i) std::cout << ','; std::cout << perGene[(size_t)i]; }
	std::cout << "]}\n";
	return 0;
}
