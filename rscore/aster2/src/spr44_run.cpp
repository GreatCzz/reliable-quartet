// spr44_run: 用精确公式驱动 spr44 搜索（每轮取"全局最优一步"）
//
// 每轮对每个候选子树 A（内部节点、非根）：
//   T2 = S 摘掉 A；备好 V_z、sumV、sumA；DFS 出前缀和 P
//   每个挂载位置 v：2*score(S') = sumA + W_w(v) + sumV + P[v]
// 取全局最大的 (A,v) 应用之。公式见 SPR44_MATH.md 第 2 节。
//
// 用法： spr44_run <species.tree> <genes.trees> <nGenes> <rounds> [threads]
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

// 返回 2*score = Σ_z W_z（不做 /2，避免整除丢半格）
long long W_of(vector<pair<int, vector<int>>> const& assign) {
	static thread_local rscore::Scratch sc;
	static thread_local vector<unsigned char> color;
	if ((int)color.size() < nTaxa) color.assign((size_t)nTaxa, 0);
	long long tot = 0;
	for (auto const& g : genes) {
		if ((int)sc.c1.size() < g.n) sc.resize(g.n);
		for (auto const& pr : assign)
			for (int t : pr.second) color[(size_t)t] = (unsigned char)pr.first;
		tot += rscore::W_colors(g, color, sc);
		for (auto const& pr : assign)
			for (int t : pr.second) color[(size_t)t] = 0;
	}
	return tot;
}

long long score2Tree(Tree const& T) {
	long long tot = 0;
	for (size_t u = 0; u < T.ch.size(); ++u) {
		if (T.taxon[u] != -1 || T.ch[u].size() != 2) continue;
		tot += W_of({{1, T.leafset[T.ch[u][0]]}, {2, T.leafset[T.ch[u][1]]}});
	}
	return tot;
}

string toNewick(Tree const& T, int u) {
	if (T.taxon[u] >= 0) return common::taxonName2ID[(size_t)T.taxon[u]];
	string r = "(";
	for (size_t i = 0; i < T.ch[u].size(); ++i) { if (i) r += ","; r += toNewick(T, T.ch[u][i]); }
	return r + ")";
}

Tree detach(Tree const& T, int A) {
	vector<char> inA(T.ch.size(), 0);
	{ vector<int> st{A}; while (!st.empty()) { int x = st.back(); st.pop_back(); inA[(size_t)x] = 1; for (int c : T.ch[x]) st.push_back(c); } }
	Tree R;
	std::function<int(int)> build = [&](int u) -> int {
		int const nu = R.addNode(T.taxon[u]);
		vector<int> kids;
		for (int c : T.ch[u]) if (!inA[(size_t)c]) kids.push_back(build(c));
		if (kids.size() == 1) { R.taxon[(size_t)nu] = -2; return kids[0]; }
		R.ch[(size_t)nu] = kids;
		for (int c : kids) R.par[(size_t)c] = nu;
		return nu;
	};
	R.root = build(T.root);
	R.par[(size_t)R.root] = -1;
	R.recomputeLeafset();
	return R;
}

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
	else { for (auto& c : S.ch[(size_t)u]) if (c == v) c = w; S.par[(size_t)w] = u; }
	S.recomputeLeafset();
	return S;
}

// 去掉死节点、重排索引
Tree compact(Tree const& T) {
	Tree R;
	std::function<int(int)> rec = [&](int u) -> int {
		int const nu = R.addNode(T.taxon[u]);
		vector<int> kids;
		for (int c : T.ch[u]) kids.push_back(rec(c));
		if (kids.size() == 1) { R.taxon[(size_t)nu] = -2; return kids[0]; }
		R.ch[(size_t)nu] = kids;
		for (int c : kids) R.par[(size_t)c] = nu;
		return nu;
	};
	R.root = rec(T.root);
	R.par[(size_t)R.root] = -1;
	R.recomputeLeafset();
	return R;
}

} // namespace

int main(int argc, char** argv) {
	if (argc < 5) { cerr << "usage: spr44_run <species.tree> <genes.trees> <nGenes> <rounds> [threads]\n"; return 1; }
	string spFile = argv[1], gFile = argv[2];
	int const nG = atoi(argv[3]);
	int const nRounds = atoi(argv[4]);
	int const nThreads = argc > 5 ? atoi(argv[5]) : 8;

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
	nTaxa = (int)common::taxonName2ID.nTaxa();

	auto const t0 = chrono::steady_clock::now();
	long long cur = score2Tree(T) / 2;
	cerr << "[spr44] taxa=" << nTaxa << " genes=" << genes.size() << " threads=" << nThreads
	     << " start=" << cur << "\n";

	for (int it = 0; it < nRounds; ++it) {
		auto const r0 = chrono::steady_clock::now();
		int const nNode = (int)T.ch.size();
		long long best = cur * 2;           // 用 2*score 尺度比较
		int bestA = -1, bestV = -1;
		Tree bestT;

		#pragma omp parallel num_threads(nThreads)
		{
			long long lbest = cur * 2;
			int lA = -1, lV = -1;
			Tree lT;
			#pragma omp for schedule(dynamic, 1)
			for (int A = 0; A < nNode; ++A) {
				if (T.taxon[(size_t)A] != -1 || A == T.root || T.ch[(size_t)A].size() != 2) continue;
				Tree const T2 = detach(T, A);
				size_t const n = T2.ch.size();
				vector<int> const X = T.leafset[(size_t)A];

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
				vector<long long> V(n, 0);
				long long sumV = 0;
				for (size_t z = 0; z < n; ++z) {
					if (T2.taxon[z] != -1 || T2.ch[z].size() != 2) continue;
					V[z] = W_of({{1, T2.leafset[T2.ch[z][0]]}, {2, T2.leafset[T2.ch[z][1]]}});
					sumV += V[z];
				}
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
				long long const base = sumA + sumV;
				for (size_t v = 0; v < n; ++v) {
					if (T2.taxon[v] == -2 || (int)v == T2.root) continue;
					if (T2.par[v] < 0 || T2.ch[(size_t)T2.par[v]].size() != 2) continue;
					long long const s = base + W_of({{1, X}, {2, T2.leafset[v]}}) + P[v];
					if (s > lbest) {
						lbest = s; lA = A; lV = (int)v;
						lT = insertA(T2, (int)v, T, A);   // 保留候选树，胜出后直接用
					}
				}
			}
			#pragma omp critical
			{
				if (lbest > best && lA >= 0) { best = lbest; bestA = lA; bestV = lV; bestT = std::move(lT); }
			}
		}

		auto const r1 = chrono::steady_clock::now();
		if (bestA < 0) { cerr << "[spr44] converged at round " << it << "\n"; break; }
		T = compact(bestT);
		long long const got = score2Tree(T) / 2;
		cerr << "[spr44] round " << it << ": " << cur << " -> " << got
		     << " (+" << (got - cur) << ")  A=" << bestA
		     << "  " << chrono::duration<double>(r1 - r0).count() << " s\n";
		if (got != best / 2 && got != (best + 1) / 2)
			cerr << "[spr44]   WARNING: 公式预测 " << best << " (2x) vs 实测 " << got * 2 << "\n";
		cur = got;
	}
	auto const t1 = chrono::steady_clock::now();
	cerr << "[spr44] final = " << cur << "  time = " << chrono::duration<double>(t1 - t0).count() << " s\n";
	cout << toNewick(T, T.root) << ";\n";
	return 0;
}
