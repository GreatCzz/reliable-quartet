// spr44_bf: spr44 的暴力参照实现 —— 对每个候选剪枝子树、每个挂载位置，构造候选树并整树精确打分
// 用途：作为"每轮取全局最优一步"的正确性基准，以及后续增量版的对照。
#include "rscore.hpp"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <functional>

using namespace std;

namespace {

struct Tree {
	vector<vector<int>> ch;
	vector<int> par;
	vector<int> taxon;
	vector<vector<int>> leafset;   // 每个节点的叶子（taxon id）列表
	int root = 0;
	int nTaxa = 0;

	int addLeaf(int tx) {
		int id = (int)ch.size();
		ch.push_back({}); par.push_back(-1); taxon.push_back(tx);
		return id;
	}
	int addInternal(int a, int b) {
		int id = (int)ch.size();
		ch.push_back({a, b}); par.push_back(-1); taxon.push_back(-1);
		par[a] = id; par[b] = id;
		return id;
	}
	// 从 root DFS，重建干净数组（丢掉不可达节点、修正 par）
	void compact() {
		vector<int> order, st{root};
		while (!st.empty()) { int u = st.back(); st.pop_back(); order.push_back(u);
			for (int c : ch[u]) if (taxon[c] != -2) st.push_back(c); }
		vector<int> id(ch.size(), -1);
		for (size_t i = 0; i < order.size(); ++i) id[(size_t)order[i]] = (int)i;
		vector<vector<int>> nch(order.size());
		vector<int> npar(order.size(), -1), ntax(order.size(), -1);
		for (size_t i = 0; i < order.size(); ++i) {
			int const u = order[i];
			ntax[i] = (taxon[u] == -2 ? -1 : taxon[u]);
			for (int c : ch[u]) if (id[(size_t)c] >= 0 && taxon[c] != -2) { nch[i].push_back(id[(size_t)c]); npar[id[(size_t)c]] = (int)i; }
		}
		ch = std::move(nch); par = std::move(npar); taxon = std::move(ntax);
		root = 0;
		recomputeLeafset();
	}

	void recomputeLeafset() {
		leafset.assign(ch.size(), {});
		std::function<void(int)> dfs = [&](int u) {
			if (taxon[u] >= 0) { leafset[(size_t)u] = {taxon[u]}; return; }
			leafset[(size_t)u].clear();
			for (int c : ch[u]) { dfs(c); auto& s = leafset[(size_t)u]; s.insert(s.end(), leafset[(size_t)c].begin(), leafset[(size_t)c].end()); }
		};
		if (!ch.empty()) dfs(root);
	}
};

struct Loader {
	vector<rscore::GeneTree> genes;
	Tree T;

	int build(newick::Newick::Node const* u) {
		vector<int> kids;
		for (auto* c : u->children) kids.push_back(build(c));
		if (kids.empty()) return T.addLeaf((int)common::taxonName2ID[u->name]);
		return T.addInternal(kids[0], kids.size() > 1 ? kids[1] : kids[0]);
	}
};

string toNewick(Tree const& T, int u) {
	if (T.taxon[u] >= 0) return common::taxonName2ID[(size_t)T.taxon[u]];
	string r = "(";
	for (size_t i = 0; i < T.ch[u].size(); ++i) { if (i) r += ","; r += toNewick(T, T.ch[u][i]); }
	return r + ")";
}

// 整树精确打分（用 rscore::W_colors，与 solver 逐值一致）
long long scoreTree(Tree const& T, vector<rscore::GeneTree> const& genes) {
	long long total = 0;
	for (auto const& g : genes) {
		rscore::Scratch sc; sc.resize(g.n);
		vector<unsigned char> color(common::taxonName2ID.nTaxa(), 0);
		for (size_t u = 0; u < T.ch.size(); ++u) {
			if (T.taxon[u] >= 0 || T.ch[u].size() != 2 || T.taxon[u] == -2) continue;
			for (int t : T.leafset[T.ch[u][0]]) color[(size_t)t] = 1;
			for (int t : T.leafset[T.ch[u][1]]) color[(size_t)t] = 2;
			total += rscore::W_colors(g, color, sc);
			for (int t : T.leafset[T.ch[u][0]]) color[(size_t)t] = 0;
			for (int t : T.leafset[T.ch[u][1]]) color[(size_t)t] = 0;
		}
	}
	return total / 2;
}

} // namespace

int main(int argc, char** argv) {
	string spFile = argv[1], gFile = argv[2];
	int nG = argc > 3 ? atoi(argv[3]) : 20;
	int maxIter = argc > 4 ? atoi(argv[4]) : 20;
	newick::Newick sp;
	{ ifstream fin(spFile); string line; while (getline(fin, line)) if (line.find_first_not_of(" \t\r\n") != string::npos) { sp = newick::Newick(line); break; } }
	Loader L;
	L.T.root = L.build(sp.root());
	L.T.recomputeLeafset();
	L.T.nTaxa = (int)common::taxonName2ID.nTaxa();
	{ ifstream fin(gFile); string line; int cnt = 0;
	  while (getline(fin, line) && cnt < nG) {
		if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
		newick::Newick nw(line); L.genes.push_back(rscore::buildGeneTree(nw.root())); ++cnt; } }

	auto t0 = chrono::steady_clock::now();
	long long cur = scoreTree(L.T, L.genes);
	cerr << "[spr44bf] taxa=" << L.T.nTaxa << " genes=" << L.genes.size() << " start=" << cur << "\n";
	int nNode = (int)L.T.ch.size();
	for (int it = 0; it < maxIter; ++it) {
		long long best = cur; int bestA = -1, bestV = -1; Tree bestT2;
		// 对每个候选剪枝子树 A（内部节点、非根）
		for (int A = 0; A < nNode; ++A) {
			if (L.T.taxon[A] != -1 || A == L.T.root) continue;   // 只取活着的内部节点
			// 对每个挂载点 v（不在 subtree(A) 内）
			for (int v = 0; v < nNode; ++v) {
				if (v == A || v == L.T.root || L.T.taxon[v] == -2) continue;
				bool inside = false;
				for (int t : L.T.leafset[A]) if (find(L.T.leafset[v].begin(), L.T.leafset[v].end(), t) != L.T.leafset[v].end()) { inside = true; break; }
				if (inside) continue;
				// 构造：摘 A、在 v 上方插 A
				Tree T2 = L.T;
				int p = T2.par[A], g = T2.par[p];
				int sib = (T2.ch[p][0] == A) ? T2.ch[p][1] : T2.ch[p][0];
				if (g < 0) { T2.root = sib; T2.par[sib] = -1; }
				else { for (auto& c : T2.ch[g]) if (c == p) c = sib; T2.par[sib] = g; }
				int u = T2.par[v];
				int w = (int)T2.ch.size();
				T2.ch.push_back({A, v}); T2.par.push_back(u); T2.taxon.push_back(-1);
				if (u < 0) { T2.root = w; T2.par[w] = -1; }
				else { for (auto& c : T2.ch[u]) if (c == v) c = w; T2.par[w] = u; }
				T2.par[A] = w; T2.par[v] = w;
				T2.taxon[p] = -2;   // 旧的 p 已被抑制，必须标死（否则 A 会被重复计数）
				T2.recomputeLeafset();
				long long s = scoreTree(T2, L.genes);
				if (s > best) { best = s; bestA = A; bestV = v; bestT2 = std::move(T2); }
			}
		}
		if (bestA < 0) { cerr << "[spr44bf] converged at iter " << it << "\n"; break; }
		// 应用：直接采用胜出的 T2（避免重复实现移动逻辑）
		L.T = std::move(bestT2);
		L.T.compact();
		{
			string const nw = toNewick(L.T, L.T.root) + ";";
			cerr << "[spr44bf] DBG emitted=" << nw.substr(0, 120) << "\n";
			cerr << "[spr44bf] DBG rescore(in-memory)=" << scoreTree(L.T, L.genes) << "\n";
		}
		cerr << "[spr44bf] iter " << it << ": " << cur << " -> " << best << " (+" << (best - cur) << ")\n";
		cur = best;
	}
	auto t1 = chrono::steady_clock::now();
	cerr << "[spr44bf] final = " << cur << "  time = " << chrono::duration<double>(t1 - t0).count() << " s\n";
	cout << toNewick(L.T, L.T.root) << ";\n";
	return 0;
}
