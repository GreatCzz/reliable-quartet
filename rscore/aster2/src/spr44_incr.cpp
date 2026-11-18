// spr44_incr_check: 用"增量公式"算一个候选子树 A 在每个挂载位置的 Δ_new，并与暴力（构造+整树打分）对拍
//   Δ_new(e) = [ W(A, subtree(v), rest) + P(u) ] / 2,   P(u) = Σ_{z on root->u} L(z, child)
//   L(z,c)   = W(A ∪ subtree(c), sibling(c), rest) − W(subtree(c), sibling(c), rest)
#include "rscore.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <functional>
#include <vector>
using namespace std;

namespace {
struct Tree {
	vector<vector<int>> ch; vector<int> par, taxon; vector<vector<int>> leafset;
	int root = 0;
	void recomputeLeafset() {
		leafset.assign(ch.size(), {});
		std::function<void(int)> dfs = [&](int u) {
			if (taxon[u] >= 0) { leafset[(size_t)u] = {taxon[u]}; return; }
			leafset[(size_t)u].clear();
			for (int c : ch[u]) { dfs(c); auto& s = leafset[(size_t)u]; s.insert(s.end(), leafset[(size_t)c].begin(), leafset[(size_t)c].end()); }
		};
		if (!ch.empty()) dfs(root);
	}
	void compact() {
		vector<int> order, st{root};
		while (!st.empty()) { int u = st.back(); st.pop_back(); order.push_back(u); for (int c : ch[u]) if (taxon[c] != -2) st.push_back(c); }
		vector<int> id(ch.size(), -1);
		for (size_t i = 0; i < order.size(); ++i) id[(size_t)order[i]] = (int)i;
		vector<vector<int>> nch(order.size()); vector<int> npar(order.size(), -1), ntax(order.size(), -1);
		for (size_t i = 0; i < order.size(); ++i) { int u = order[i]; ntax[i] = (taxon[u] == -2 ? -1 : taxon[u]);
			for (int c : ch[u]) if (id[(size_t)c] >= 0 && taxon[c] != -2) { nch[i].push_back(id[(size_t)c]); npar[id[(size_t)c]] = (int)i; } }
		ch = std::move(nch); par = std::move(npar); taxon = std::move(ntax); root = 0; recomputeLeafset();
	}
};

vector<rscore::GeneTree> genes;
int nTaxa = 0;

long long W_of(vector<pair<int, vector<int>>> const& assign) {   // (mask, leafset) 列表
	long long tot = 0;
	for (auto const& g : genes) {
		rscore::Scratch sc; sc.resize(g.n);
		vector<unsigned char> color((size_t)nTaxa, 0);
		for (auto const& pr : assign) for (int t : pr.second) color[(size_t)t] = (unsigned char)pr.first;
		tot += rscore::W_colors(g, color, sc);
	}
	return tot;
}

long long scoreTree(Tree const& T) {
	long long tot = 0;
	for (size_t u = 0; u < T.ch.size(); ++u) {
		if (T.taxon[u] >= 0 || T.ch[u].size() != 2) continue;
		tot += W_of({{1, T.leafset[T.ch[u][0]]}, {2, T.leafset[T.ch[u][1]]}});
	}
	return tot / 2;
}

string toNewick(Tree const& T, int u) {
	if (T.taxon[u] >= 0) return common::taxonName2ID[(size_t)T.taxon[u]];
	string r = "(";
	for (size_t i = 0; i < T.ch[u].size(); ++i) { if (i) r += ","; r += toNewick(T, T.ch[u][i]); }
	return r + ")";
}
}


// 对给定 A，返回每个合法位置 v 的相对 Δ（相对 A 的原位置），并返回最优
struct PlaceRes { long long bestDelta = 0; int bestV = -1; };

PlaceRes placeOfA(Tree const& T, int A) {
	Tree Tp = T;
	int p = Tp.par[A], g = Tp.par[p];
	int sib = (Tp.ch[p][0] == A) ? Tp.ch[p][1] : Tp.ch[p][0];
	if (g < 0) { Tp.root = sib; Tp.par[sib] = -1; }
	else { for (auto& c : Tp.ch[g]) if (c == p) c = sib; Tp.par[sib] = g; }
	Tp.taxon[p] = -2; Tp.ch[p].clear(); Tp.recomputeLeafset();
	vector<int> Ax = T.leafset[(size_t)A];
	size_t const n = Tp.ch.size();
	vector<long long> Worig(n, 0), P(n, 0);
	for (size_t z = 0; z < n; ++z) {
		if (Tp.taxon[z] != -1 || Tp.ch[z].size() != 2) continue;
		Worig[z] = W_of({{1, Tp.leafset[Tp.ch[z][0]]}, {2, Tp.leafset[Tp.ch[z][1]]}});
	}
	std::function<void(int)> dfs = [&](int z) {
		for (int c : Tp.ch[z]) {
			int other = (Tp.ch[z][0] == c) ? Tp.ch[z][1] : Tp.ch[z][0];
			vector<int> merged = Tp.leafset[(size_t)c];
			merged.insert(merged.end(), Ax.begin(), Ax.end());
			P[(size_t)c] = P[(size_t)z] + (W_of({{2, merged}, {1, Tp.leafset[(size_t)other]}}) - Worig[z]);
			dfs(c);
		}
	};
	dfs(Tp.root);
	// 原位置（挂到 sib 上方）的值
	long long const atOld = (W_of({{1, Ax}, {2, Tp.leafset[(size_t)sib]}}) + P[(size_t)sib]) / 2;
	PlaceRes res;
	for (size_t v = 0; v < n; ++v) {
		if ((int)v == Tp.root || Tp.taxon[v] == -2 || Tp.leafset[v].empty()) continue;
		long long const val = (W_of({{1, Ax}, {2, Tp.leafset[v]}}) + P[(size_t)v]) / 2;
		if (val - atOld > res.bestDelta) { res.bestDelta = val - atOld; res.bestV = (int)v; }
	}
	return res;
}

int main(int argc, char** argv) {
	string spFile = argv[1], gFile = argv[2];
	int nG = argc > 3 ? atoi(argv[3]) : 20;
	int maxIter = argc > 4 ? atoi(argv[4]) : 10;
	newick::Newick sp;
	{ ifstream fin(spFile); string line; while (getline(fin, line)) if (line.find_first_not_of(" \t\r\n") != string::npos) { sp = newick::Newick(line); break; } }
	Tree T;
	{ std::function<int(newick::Newick::Node const*)> build = [&](newick::Newick::Node const* u) -> int {
		vector<int> kids; for (auto* c : u->children) kids.push_back(build(c));
		int id = (int)T.ch.size(); T.ch.push_back({}); T.par.push_back(-1); T.taxon.push_back(-1);
		if (kids.empty()) T.taxon[id] = (int)common::taxonName2ID[u->name];
		else { T.ch[id] = kids; for (int c : kids) T.par[c] = id; }
		return id; };
	  T.root = build(sp.root()); }
	T.recomputeLeafset(); nTaxa = (int)common::taxonName2ID.nTaxa();
	{ ifstream fin(gFile); string line; int c = 0;
	  while (getline(fin, line) && c < nG) { if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
		newick::Newick nw(line); genes.push_back(rscore::buildGeneTree(nw.root())); ++c; } }

	auto t0 = chrono::steady_clock::now();
	long long cur = scoreTree(T);
	cerr << "[incr] taxa=" << nTaxa << " genes=" << genes.size() << " start=" << cur << "\n";
	for (int it = 0; it < maxIter; ++it) {
		long long bestD = 0; int bestA = -1, bestV = -1;
		for (size_t A = 0; A < T.ch.size(); ++A) {
			if (T.taxon[A] != -1 || (int)A == T.root || T.ch[A].size() != 2) continue;
			PlaceRes const r = placeOfA(T, (int)A);
			if (r.bestDelta > bestD) { bestD = r.bestDelta; bestA = (int)A; bestV = r.bestV; }
		}
		if (bestA < 0) { cerr << "[incr] converged at iter " << it << "\n"; break; }
		// 应用：摘 bestA，挂到 bestV 上方
		int p = T.par[bestA], g = T.par[p];
		int sib = (T.ch[p][0] == bestA) ? T.ch[p][1] : T.ch[p][0];
		if (g < 0) { T.root = sib; T.par[sib] = -1; }
		else { for (auto& c : T.ch[g]) if (c == p) c = sib; T.par[sib] = g; }
		int u = T.par[bestV];
		int w = (int)T.ch.size();
		T.ch.push_back({bestA, bestV}); T.par.push_back(u); T.taxon.push_back(-1);
		if (u < 0) { T.root = w; T.par[w] = -1; }
		else { for (auto& c : T.ch[u]) if (c == bestV) c = w; T.par[w] = u; }
		T.par[bestA] = w; T.par[bestV] = w;
		T.taxon[p] = -2; T.ch[p].clear();
		T.recomputeLeafset(); T.compact();
		cur += bestD;
		cerr << "[incr] iter " << it << ":  A=" << bestA << " v=" << bestV << " +" << bestD << " -> " << cur << "\n";
	}
	auto t1 = chrono::steady_clock::now();
	long long const fin = scoreTree(T);
	cerr << "[incr] final(预测)=" << cur << "  final(重算)=" << fin << "  time=" << chrono::duration<double>(t1 - t0).count() << " s\n";
	cout << toNewick(T, T.root) << ";\n";
	return 0;
}
