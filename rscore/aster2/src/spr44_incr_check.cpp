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
	vector<char> reach(T.ch.size(), 0);
	{ vector<int> st{T.root}; while (!st.empty()) { int u = st.back(); st.pop_back(); reach[(size_t)u] = 1; for (int c : T.ch[u]) st.push_back(c); } }
	for (size_t u = 0; u < T.ch.size(); ++u) {
		if (!reach[u]) continue;
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

int main(int argc, char** argv) {
	string spFile = argv[1], gFile = argv[2];
	int nG = argc > 3 ? atoi(argv[3]) : 20;
	newick::Newick sp;
	{ ifstream fin(spFile); string line; while (getline(fin, line)) if (line.find_first_not_of(" \t\r\n") != string::npos) { sp = newick::Newick(line); break; } }
	Tree T; vector<int> taxonOf;
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

	// 选一个候选 A：第一个非根内部节点
	int A = -1;
	if (argc > 4) A = atoi(argv[4]);
	else for (size_t u = 0; u < T.ch.size(); ++u) if (T.taxon[u] < 0 && (int)u != T.root && T.ch[u].size() == 2) { A = (int)u; break; }
	if (A < 0 || T.taxon[A] >= 0 || (int)A == T.root) { cerr << "no candidate\n"; return 1; }

	// T' = S - A
	Tree Tp = T;
	{ int p = Tp.par[A], g = Tp.par[p];
	  int sib = (Tp.ch[p][0] == A) ? Tp.ch[p][1] : Tp.ch[p][0];
	  if (g < 0) { Tp.root = sib; Tp.par[sib] = -1; } else { for (auto& c : Tp.ch[g]) if (c == p) c = sib; Tp.par[sib] = g; }
	  Tp.taxon[p] = -2; Tp.ch[p].clear(); Tp.recomputeLeafset(); }
	long long const base = scoreTree(Tp);
	vector<int> Ax = T.leafset[A];

	// W_orig(z) 与 L(z,c)
	size_t const n = Tp.ch.size();
	vector<long long> Worig(n, 0), P(n, 0);
	{
		vector<char> reach(n, 0);
		vector<int> st{Tp.root};
		while (!st.empty()) { int u = st.back(); st.pop_back(); reach[(size_t)u] = 1; for (int c : Tp.ch[u]) st.push_back(c); }
		for (size_t z = 0; z < n; ++z) {
			if (!reach[z]) continue;
			if (Tp.taxon[z] != -1 || Tp.ch[z].size() != 2) continue;
			Worig[z] = W_of({{1, Tp.leafset[Tp.ch[z][0]]}, {2, Tp.leafset[Tp.ch[z][1]]}});
		}
	}
	// 前缀和（root 为 0）
	std::function<void(int)> dfs = [&](int z) {
		for (int c : Tp.ch[z]) {
			int other = (Tp.ch[z][0] == c) ? Tp.ch[z][1] : Tp.ch[z][0];
			vector<int> merged = Tp.leafset[(size_t)c];
			merged.insert(merged.end(), Ax.begin(), Ax.end());
			long long const Wa = W_of({{2, merged}, {1, Tp.leafset[(size_t)other]}});
			P[(size_t)c] = P[(size_t)z] + (Wa - Worig[z]);
			dfs(c);
		}
	};
	dfs(Tp.root);

	// 对比：每个位置 v 的 Δ_new
	int bad = 0, tot = 0; long long pred0 = 0, got0 = 0; bool first = true;
	for (size_t v = 0; v < n; ++v) {
		if ((int)v == Tp.root || Tp.taxon[v] == -2 || Tp.leafset[v].empty()) continue;   // 跳过死节点
		int u = Tp.par[v];
		long long const Ww = W_of({{1, Ax}, {2, Tp.leafset[v]}});
		long long const pred = (Ww + P[(size_t)v]) / 2;     // 路径和应含 u 自身的变化 ⇒ 用 P[v]
		// 暴力：把 A 挂到 v 上方，整树打分
		Tree T2 = Tp;
		int w = (int)T2.ch.size();
		T2.ch.push_back({A, (int)v}); T2.par.push_back(u); T2.taxon.push_back(-1);
		if (u < 0) { T2.root = w; T2.par[w] = -1; } else { for (auto& c : T2.ch[u]) if (c == (int)v) c = w; T2.par[w] = u; }
		T2.par[A] = w; T2.par[(size_t)v] = w;
		T2.recomputeLeafset();
		long long const got = scoreTree(T2) - base;
		++tot;
		if (first) { pred0 = pred; got0 = got; first = false; }   // 只比较相对值（A 内部节点的分数是放置无关常数）
		if (pred - pred0 != got - got0) { ++bad; if (bad <= 8) cerr << "  MISMATCH v=" << v << " sz(v)=" << Tp.leafset[v].size() << " pred=" << (pred - pred0) << " got=" << (got - got0) << "\n"; }
	}
	cerr << "[incr-check] A=" << A << " leaves=" << Ax.size() << " positions=" << tot << " mismatches=" << bad << "\n";
	return bad == 0 ? 0 : 1;
}
