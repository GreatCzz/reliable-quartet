// spr44_formula: 验证 spr44 的精确增量公式（两阶段分解，无残差）
//
// 记号：S 原树，A 待移动子树（叶子集 X），p = par(A)，g = par(p)，sib = p 的另一个孩子
//       T2 = S − A（p 摘掉、g 抑制），score(T) = Σ_z W_z / 2
//
// 分解（引理：detach 与 insert 是两次独立的改动，各自树内成立）：
//   Δ(A,e) = Δ_old(A) + Δ_new(A,e)
//   Δ_old(A)      = score(T2) − score(S)                      （与 e 无关）
//   Δ_new(A,e)    = [ W(X, sub(v), rest_w) + P(u) ] / 2
//     P(root) = 0
//     P(c)    = P(z) + L(z,c)     对 T2 的每条有向边 (z → c)
//     L(z,c)  = W(sub(c) ∪ X, sub(sib_z), rest_z) − W(sub(c), sub(sib_z), rest_z)
//     rest_w  = T2 中 w 之外的全部类群 = rest_z(u) ∪ sub(sib_z(u))
//
// 关键点：L 与 W_w 里用到的 "rest" 全部取自 **T2**（A 已摘掉的那棵树），
//         不能一个用 S、一个用 T2 —— 之前 11~38 的残差就来自这里。
#include "rscore.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

namespace {

struct Tree {
	vector<vector<int>> ch;
	vector<int> par;
	vector<int> taxon;          // >=0 叶子；-1 内部；-2 死节点
	vector<vector<int>> leafset;
	int root = 0;

	int addNode(int tx) {
		int id = (int)ch.size();
		ch.push_back({}); par.push_back(-1); taxon.push_back(tx);
		return id;
	}
	void recomputeLeafset() {
		leafset.assign(ch.size(), {});
		std::function<void(int)> dfs = [&](int u) {
			if (taxon[u] >= 0) { leafset[(size_t)u] = {taxon[u]}; return; }
			leafset[(size_t)u].clear();
			for (int c : ch[u]) {
				dfs(c);
				auto& s = leafset[(size_t)u];
				s.insert(s.end(), leafset[(size_t)c].begin(), leafset[(size_t)c].end());
			}
		};
		if (!ch.empty()) dfs(root);
	}
};

vector<rscore::GeneTree> genes;
int nTaxa = 0;

// hide：这些类群在 gene tree 里"当作不存在"（颜色 3 不落在任何 mask 里）。
// 只有 T2 上取 W 时才需要 hide = X —— 否则 X 会被算进第三组 "rest"，得到 S 上的值而不是 T2 上的值。
long long W_of(vector<pair<int, vector<int>>> const& assign, vector<int> const* hide = nullptr) {
	long long tot = 0;
	for (auto const& g : genes) {
		rscore::Scratch sc; sc.resize(g.n);
		vector<unsigned char> color((size_t)nTaxa, 0);
		if (hide) for (int t : *hide) color[(size_t)t] = 3;
		for (auto const& pr : assign) for (int t : pr.second) color[(size_t)t] = (unsigned char)pr.first;
		tot += rscore::W_colors(g, color, sc);
	}
	return tot;
}

// 注意：返回 2*score = Σ_z W_z，**不做 /2**。除不尽时会丢半个单位，
// 而一切比较都应该在 "2*score" 这个整数尺度上做，才不会引入 ±1/±0.5 的伪残差。
long long score2Tree(Tree const& T, vector<int> const* hide = nullptr) {
	long long tot = 0;
	for (size_t u = 0; u < T.ch.size(); ++u) {
		if (T.taxon[u] != -1 || T.ch[u].size() != 2) continue;
		tot += W_of({{1, T.leafset[T.ch[u][0]]}, {2, T.leafset[T.ch[u][1]]}}, hide);
	}
	return tot;
}

string toNewick(Tree const& T, int u) {
	if (T.taxon[u] >= 0) return common::taxonName2ID[(size_t)T.taxon[u]];
	string r = "(";
	for (size_t i = 0; i < T.ch[u].size(); ++i) { if (i) r += ","; r += toNewick(T, T.ch[u][i]); }
	return r + ")";
}

// 摘掉 A（含抑制单叉节点），返回全新的树
Tree detach(Tree const& T, int A) {
	vector<char> inA(T.ch.size(), 0);
	{ vector<int> st{A}; while (!st.empty()) { int x = st.back(); st.pop_back(); inA[(size_t)x] = 1; for (int c : T.ch[x]) st.push_back(c); } }
	Tree R;
	std::function<int(int)> build = [&](int u) -> int {
		int const nu = R.addNode(T.taxon[u]);
		vector<int> kids;
		for (int c : T.ch[u]) if (!inA[(size_t)c]) kids.push_back(build(c));
		if (kids.size() == 1) {              // 抑制：本节点标死，把孩子直接顶上去
			R.taxon[(size_t)nu] = -2;
			return kids[0];
		}
		R.ch[(size_t)nu] = kids;
		for (int c : kids) R.par[(size_t)c] = nu;
		return nu;
	};
	R.root = build(T.root);
	R.par[(size_t)R.root] = -1;
	R.recomputeLeafset();
	return R;
}

// 把 A 的子树原样复制，插到 T2 的边 (u=par(v), v) 上
Tree insertA(Tree const& T2, int v, Tree const& T, int A) {
	Tree S = T2;
	std::function<int(int)> copySub = [&](int u) -> int {
		int const nu = S.addNode(T.taxon[u]);
		if (T.taxon[u] < 0) {
			vector<int> kids;
			for (int c : T.ch[u]) kids.push_back(copySub(c));
			S.ch[(size_t)nu] = kids;
			for (int c : kids) S.par[(size_t)c] = nu;
		}
		return nu;
	};
	int const a2 = copySub(A);
	int const w = S.addNode(-1);
	int const u = S.par[(size_t)v];
	S.ch[(size_t)w] = {a2, v};
	S.par[(size_t)a2] = w; S.par[(size_t)v] = w;
	if (u < 0) { S.root = w; S.par[(size_t)w] = -1; }
	else {
		for (auto& c : S.ch[(size_t)u]) if (c == v) c = w;
		S.par[(size_t)w] = u;
	}
	S.recomputeLeafset();
	return S;
}

} // namespace

int main(int argc, char** argv) {
	string spFile = argv[1], gFile = argv[2];
	int nG = argc > 3 ? atoi(argv[3]) : 20;
	newick::Newick sp;
	{ ifstream fin(spFile); string line; while (getline(fin, line)) if (line.find_first_not_of(" \t\r\n") != string::npos) { sp = newick::Newick(line); break; } }
	Tree T;
	{
		std::function<int(newick::Newick::Node const*)> build = [&](newick::Newick::Node const* u) -> int {
			vector<int> kids;
			for (auto* c : u->children) kids.push_back(build(c));
			int const id = (int)T.ch.size();
			T.ch.push_back({}); T.par.push_back(-1); T.taxon.push_back(-1);
			if (kids.empty()) T.taxon[(size_t)id] = (int)common::taxonName2ID[u->name];
			else { T.ch[(size_t)id] = kids; for (int c : kids) T.par[(size_t)c] = id; }
			return id;
		};
		T.root = build(sp.root());
	}
	T.recomputeLeafset();
	{ ifstream fin(gFile); string line; int c = 0;
	  while (getline(fin, line) && c < nG) {
		if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
		newick::Newick nw(line); genes.push_back(rscore::buildGeneTree(nw.root())); ++c; } }
	nTaxa = (int)common::taxonName2ID.nTaxa();   // 读完所有基因树后再定颜色数组大小

	auto const t0 = chrono::steady_clock::now();
	long long const scoreS = score2Tree(T);

	int nCand = 0, nPos = 0, nBad = 0;
	long long maxAbs = 0;
	for (int A = 0; A < (int)T.ch.size(); ++A) {
		if (T.taxon[(size_t)A] != -1 || A == T.root || T.ch[(size_t)A].size() != 2) continue;
		Tree const T2 = detach(T, A);


		vector<int> const X = T.leafset[(size_t)A];
		size_t const n = T2.ch.size();

		// Worig[z] = W_z(T2)（X 从基因树里消失）；V[z] = W_z(T2) 但 X 留在 rest 里
		vector<long long> V(n, 0);
		for (size_t z = 0; z < n; ++z) {
			if (T2.taxon[z] != -1 || T2.ch[z].size() != 2) continue;
			auto const& a = T2.leafset[T2.ch[z][0]];
			auto const& b = T2.leafset[T2.ch[z][1]];
			V[z] = W_of({{1, a}, {2, b}});
		}

		// 前缀和 P'(c) = P'(z) + (Wa(z,c) − V[z])
		vector<long long> P(n, 0);
		{
			std::function<void(int)> dfs = [&](int z) {
				for (int c : T2.ch[(size_t)z]) {
					int const other = (T2.ch[(size_t)z][0] == c) ? T2.ch[(size_t)z][1] : T2.ch[(size_t)z][0];
					vector<int> merged = T2.leafset[(size_t)c];
					merged.insert(merged.end(), X.begin(), X.end());
					long long const Wa = W_of({{2, merged}, {1, T2.leafset[(size_t)other]}});
					P[(size_t)c] = P[(size_t)z] + (Wa - V[(size_t)z]);
					dfs(c);
				}
			};
			dfs(T2.root);
		}
		long long const sumV = [&] { long long t = 0; for (size_t z = 0; z < n; ++z) t += V[z]; return t; }();

		// A 自身内部节点：它们在 T2 里不存在，但在 S' 里回来了，
		// 且其三分划与挂在哪里无关（rest 恒为 A 之外的类群）→ 整块是 A 的常数
		long long sumA = 0;
		{
			vector<int> st{A};
			while (!st.empty()) {
				int const z = st.back(); st.pop_back();
				if (T.taxon[(size_t)z] == -1 && T.ch[(size_t)z].size() == 2)
					sumA += W_of({{1, T.leafset[T.ch[(size_t)z][0]]}, {2, T.leafset[T.ch[(size_t)z][1]]}});
				for (int c : T.ch[(size_t)z]) st.push_back(c);
			}
		}

		++nCand;
		for (size_t v = 0; v < n; ++v) {
			if (T2.taxon[v] == -2) continue;
			if ((int)v == T2.root) continue;                      // 根上方挂载：另算，这里先跳过
			int const u = T2.par[v];
			if (u < 0 || T2.ch[(size_t)u].size() != 2) continue;
			// rest_w = T2 中 (sub(v) ∪ X) 之外的全部类群 = R_u ∪ sub(sib(u))；
			// W_colors 的第三组就是"其余"，故只需给前两组即可
			long long const Ww = W_of({{1, X}, {2, T2.leafset[v]}});
			long long const pred2 = -scoreS + sumA + Ww + sumV + P[v];
			++nPos;
			if (getenv("DUMP")) {
				static string xkey; xkey.clear();
				{ vector<int> t = X; sort(t.begin(), t.end()); for (int z : t) { if (!xkey.empty()) xkey += ","; xkey += common::taxonName2ID[(size_t)z]; } }
				cout << xkey << "\t" << pred2 << "\n";
			}
			if (getenv("SKIPBRUTE")) { continue; }
			Tree const S2 = insertA(T2, (int)v, T, A);
			long long const got2 = score2Tree(S2) - scoreS;
			if (getenv("DBG") && nBad < 3) {
				long long lc = 0, lcT2 = 0;
				for (size_t z = 0; z < S2.ch.size(); ++z) if (S2.taxon[z] >= 0) ++lc;
				for (size_t z = 0; z < T2.ch.size(); ++z) if (T2.taxon[z] >= 0) ++lcT2;
				cerr << "  [dbg] A=" << A << " v=" << (int)v << " |X|=" << X.size()
				     << " leaves(T)=" << nTaxa << " leaves(T2)=" << lcT2 << " leaves(S2)=" << lc
				     << " score2(S2)=" << score2Tree(S2)
				     << " Ww=" << Ww << " P[v]=" << P[v] << "\n";
			}
			if (pred2 != got2) {
				++nBad;
				long long const d = pred2 - got2;
				maxAbs = max(maxAbs, d < 0 ? -d : d);
				if (nBad <= 20)
					cerr << "  MISMATCH A=" << A << " v=" << (int)v << " |X|=" << X.size()
					     << " |sub(v)|=" << T2.leafset[v].size()
					     << " pred=" << pred2 << " got=" << got2 << " diff=" << d << "\n";
			}
		}
	}
	auto const t1 = chrono::steady_clock::now();
	cerr << "[spr44-formula] taxa=" << nTaxa << " genes=" << genes.size()
	     << " cand=" << nCand << " positions=" << nPos << " mismatches=" << nBad
	     << " maxAbsDiff=" << maxAbs
	     << " scoreS=" << scoreS
	     << "  time=" << chrono::duration<double>(t1 - t0).count() << " s\n";
	return nBad == 0 ? 0 : 1;
}
