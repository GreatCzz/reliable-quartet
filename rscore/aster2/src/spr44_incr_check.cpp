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
	vector<long long> Worig(n, 0), P(n, 0), Lv(n, 0);
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
			Lv[(size_t)c] = Wa - Worig[z];
			P[(size_t)c] = P[(size_t)z] + Lv[(size_t)c];
			dfs(c);
		}
	};
	dfs(Tp.root);

	// Tp 的可达集合（位置只能取自可达节点，且不能是 A 自身）
	vector<char> reachP(n, 0);
	{ vector<int> st{Tp.root}; while (!st.empty()) { int u = st.back(); st.pop_back(); reachP[(size_t)u] = 1; for (int c : Tp.ch[u]) st.push_back(c); } }

	// 对比：每个位置 v 的 Δ_new
	int bad = 0, tot = 0; long long pred0 = 0, got0 = 0; bool first = true;
	for (size_t v = 0; v < n; ++v) {
		if ((int)v == Tp.root || !reachP[v] || (int)v == A) continue;   // 只取 T' 中可达、且不是 A 自身的节点
		int u = Tp.par[v];
		long long Ww = 0;   // 下面在 T2 构造后用真实节点计算
		if (getenv("DIAG_V") && v == (size_t)atoi(getenv("DIAG_V"))) {
			vector<int> t = Tp.leafset[v]; std::sort(t.begin(), t.end());
			cerr << "    [at-139] v=" << v << " vleaf="; for (int x : t) cerr << " " << x;
			cerr << " | Ax="; vector<int> a = Ax; std::sort(a.begin(), a.end()); for (int x : a) cerr << " " << x;
			cerr << " | Ww=" << Ww << " | genes=" << genes.size() << " nTaxa=" << nTaxa << "\n";
		}
		// 暴力：把 A 挂到 v 上方（复用旧父 p 的槽位，节点数与 Tp 相同）
		int const pA = T.par[A];
		Tree T2 = Tp;
		{
			T2.ch[pA] = {A, (int)v}; T2.taxon[pA] = -1;
			if (u < 0) { T2.root = pA; T2.par[pA] = -1; }
			else { for (auto& c : T2.ch[u]) if (c == (int)v) c = pA; T2.par[pA] = u; }
			T2.par[A] = pA; T2.par[(size_t)v] = pA;
			T2.recomputeLeafset();
		}
		long long const WwFinal = W_of({{1, T2.leafset[T2.ch[pA][0]]}, {2, T2.leafset[T2.ch[pA][1]]}});
		long long const pred = (WwFinal + P[(size_t)v]) / 2;   // 用真实新节点的孩子集合
		long long const got = scoreTree(T2) - base;
		{
			// 判断 v 是否在 A 的旧路径上（A 的祖先链）
			bool onOld = false;
			for (int z = T.par[A]; z >= 0; z = T.par[z]) if (z == (int)v) onOld = true;
			if ((pred - pred0) != (got - got0))
				cerr << "  MISMATCH v=" << v << " onOldPath=" << onOld << " sz(v)=" << Tp.leafset[v].size()
				     << " pred=" << (pred - pred0) << " got=" << (got - got0) << "\n";
		}
		if (v == 8) {
			vector<char> r2(T2.ch.size(), 0), rp(Tp.ch.size(), 0);
			{ vector<int> st{T2.root}; while (!st.empty()) { int u = st.back(); st.pop_back(); r2[(size_t)u] = 1; for (int c : T2.ch[u]) st.push_back(c); } }
			{ vector<int> st{Tp.root}; while (!st.empty()) { int u = st.back(); st.pop_back(); rp[(size_t)u] = 1; for (int c : Tp.ch[u]) st.push_back(c); } }
			long long sum = 0;
			for (size_t z = 0; z < T2.ch.size(); ++z) {
				long long w2 = 0, wp = 0;
				if (z < Tp.ch.size() && rp[z] && Tp.taxon[z] == -1 && Tp.ch[z].size() == 2)
					wp = W_of({{1, Tp.leafset[Tp.ch[z][0]]}, {2, Tp.leafset[Tp.ch[z][1]]}});
				if (r2[z] && T2.taxon[z] == -1 && T2.ch[z].size() == 2)
					w2 = W_of({{1, T2.leafset[T2.ch[z][0]]}, {2, T2.leafset[T2.ch[z][1]]}});
				if (w2 != wp) { cerr << "    node " << z << " w2=" << w2 << " wp=" << wp << " d=" << (w2 - wp) << "\n"; sum += w2 - wp; }
			}
			cerr << "  SUM/2 = " << sum / 2 << "   formula pred = " << pred << "  (Ww=" << WwFinal << " P[v]=" << P[(size_t)v] << ")\n";
		}
		if (getenv("DIAG_V") && v == (size_t)atoi(getenv("DIAG_V"))) {
			vector<char> r2(T2.ch.size(), 0), rp(Tp.ch.size(), 0);
			{ vector<int> st{T2.root}; while (!st.empty()) { int x = st.back(); st.pop_back(); r2[(size_t)x] = 1; for (int c : T2.ch[x]) st.push_back(c); } }
			{ vector<int> st{Tp.root}; while (!st.empty()) { int x = st.back(); st.pop_back(); rp[(size_t)x] = 1; for (int c : Tp.ch[x]) st.push_back(c); } }
			// A 的子树（在 T 里）
			vector<char> inA(T.ch.size(), 0);
			{ vector<int> st{A}; while (!st.empty()) { int x = st.back(); st.pop_back(); inA[(size_t)x] = 1; for (int c : T.ch[x]) st.push_back(c); } }
			bool const dbgIsA = inA[(size_t)A];
			(void)dbgIsA;
			long long sum = 0;
			for (size_t z = 0; z < T2.ch.size(); ++z) {
				long long w2 = 0, wp = 0;
				if (z < Tp.ch.size() && rp[z] && Tp.taxon[z] == -1 && Tp.ch[z].size() == 2)
					wp = W_of({{1, Tp.leafset[Tp.ch[z][0]]}, {2, Tp.leafset[Tp.ch[z][1]]}});
				if (r2[z] && T2.taxon[z] == -1 && T2.ch[z].size() == 2)
					w2 = W_of({{1, T2.leafset[T2.ch[z][0]]}, {2, T2.leafset[T2.ch[z][1]]}});
				if (w2 != wp) {
					string tag = (z == (size_t)pA) ? "W" : (inA[z] ? "inA" : "other");
					cerr << "    z=" << z << " [" << tag << "] w2=" << w2 << " wp=" << wp << " d=" << (w2 - wp) << "\n";
					sum += w2 - wp;
				}
			}
			// 打印 w 与"负项"节点的具体三分划
			for (size_t zz : {(size_t)pA, (size_t)38, (size_t)14}) {
				if (zz >= T2.ch.size()) continue;
				auto dump = [&](char const* tag, Tree const& TT) {
					if (zz >= TT.ch.size()) { cerr << "    " << tag << " z=" << zz << " (absent)\n"; return; }
					cerr << "    " << tag << " z=" << zz << " taxon=" << TT.taxon[zz] << " ch=";
					for (int c : TT.ch[zz]) cerr << c << "(" << TT.leafset[(size_t)c].size() << ") ";
					cerr << " leaves=" << TT.leafset[zz].size() << "\n";
				};
				dump("T2:", T2); dump("Tp:", Tp);
			}
			{
				auto pr = [&](char const* tag, vector<int> const& vv) {
					vector<int> t = vv; std::sort(t.begin(), t.end());
					cerr << "    " << tag << " (" << t.size() << "):"; for (int x : t) cerr << " " << x; cerr << "\n";
				};
				pr("Ax", Ax);
				if (T2.ch.size() > 14 && T2.ch[14].size() == 2) {
					pr("T2.ch[14][0]", T2.leafset[T2.ch[14][0]]);
					pr("T2.ch[14][1]", T2.leafset[T2.ch[14][1]]);
					cerr << "    W(Ax, vleaf) run1=" << W_of({{1, Ax}, {2, Tp.leafset[(size_t)v]}})
					     << " run2=" << W_of({{1, Ax}, {2, Tp.leafset[(size_t)v]}})
					     << " run3=" << W_of({{1, Ax}, {2, Tp.leafset[(size_t)v]}}) << "\n";
					cerr << "    W(ch0,ch1) =" << W_of({{1, T2.leafset[T2.ch[14][0]]}, {2, T2.leafset[T2.ch[14][1]]}}) << "\n";
				}
			}
			cerr << "  SUM/2=" << sum / 2 << "  pred=" << pred << "  Ww=" << WwFinal << "  P[v]=" << P[(size_t)v] << "\n";
			cerr << "  P-chain (root→v):";
			{ vector<int> chain; for (int z = (int)v; z >= 0; z = Tp.par[z]) chain.push_back(z); std::reverse(chain.begin(), chain.end());
			  for (int z : chain) cerr << " " << z << "(L=" << Lv[(size_t)z] << ")"; }
			cerr << "\n";
			exit(0);
		}
		++tot;
		if (first) { pred0 = pred; got0 = got; first = false; }
		bool inSib = false;
		{ int const pp = T.par[A]; int const ss = (T.ch[pp][0] == A) ? T.ch[pp][1] : T.ch[pp][0];
		  vector<char> mark(T.ch.size(), 0);
		  { vector<int> st{ss}; while (!st.empty()) { int x = st.back(); st.pop_back(); mark[(size_t)x] = 1; for (int c : T.ch[x]) st.push_back(c); } }
		  if (v < mark.size() && mark[v]) inSib = true; }
		bool ancOfP = false, pAncOfV = false;
		{ int const pp = T.par[A];
		  for (int z = pp; z >= 0; z = T.par[z]) if (z == (int)v) ancOfP = true;
		  for (int z = (int)v; z >= 0; z = T.par[z]) if (z == pp) pAncOfV = true; }
		if ((pred - pred0) != (got - got0)) cerr << "    v=" << v << " ancOfP=" << ancOfP << " pAncOfV=" << pAncOfV << " sz=" << Tp.leafset[v].size()
			<< " pred-rel=" << (pred - pred0) << " got-rel=" << (got - got0)
			<< " diff=" << ((pred - pred0) - (got - got0)) << "\n";
		if (pred - pred0 != got - got0) ++bad;
	}
	cerr << "[incr-check] A=" << A << " leaves=" << Ax.size() << " positions=" << tot << " mismatches=" << bad << "\n";
	return bad == 0 ? 0 : 1;
}
