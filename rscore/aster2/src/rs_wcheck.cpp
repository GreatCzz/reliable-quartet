// rs_wcheck: 用 rscore::W_colors 逐节点算三分划分数，与整树打分对拍
#include "rscore.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

namespace {
struct SNode { vector<int> ch; int taxon = -1; };
struct STree { vector<SNode> nd; int root = 0; vector<int> taxon; };

int build(newick::Newick::Node const* u, STree& t) {
	vector<int> kids;
	for (auto* c : u->children) kids.push_back(build(c, t));   // 先建孩子，避免 realloc 后引用失效
	int const id = (int)t.nd.size();
	t.nd.push_back({});
	t.nd[id].ch = kids;
	if (u->children.empty()) {
		t.nd[id].taxon = (int)common::taxonName2ID[u->name];
		t.taxon.push_back(t.nd[id].taxon);
	}
	return id;
}
}

int main(int argc, char** argv) {
	string spFile = argv[1], gFile = argv[2];
	int nG = argc > 3 ? atoi(argv[3]) : 100;
	newick::Newick sp(spFile == "-" ? string{} : string{});
	{ ifstream fin(spFile); string line; while (getline(fin, line)) if (line.find_first_not_of(" \t\r\n") != string::npos) { sp = newick::Newick(line); break; } }
	STree T; T.root = build(sp.root(), T);
	int n = (int)T.nd.size();

	// 子树叶子（用 taxon 位）
	vector<vector<int>> leaf(n);
	for (int u = 0; u < n; ++u) {   // 后序构造：孩子 id 更小
		if (T.nd[u].taxon >= 0) leaf[u] = {T.nd[u].taxon};
		else { for (int c : T.nd[u].ch) { leaf[u].insert(leaf[u].end(), leaf[c].begin(), leaf[c].end()); } }
	}

	vector<rscore::GeneTree> genes;
	{ ifstream fin(gFile); string line; int cnt = 0;
	  while (getline(fin, line) && cnt < nG) {
		if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
		newick::Newick nw(line); genes.push_back(rscore::buildGeneTree(nw.root())); ++cnt; } }

	long long total = 0;
	for (auto const& g : genes) {
		rscore::Scratch sc; sc.resize(g.n);
		vector<unsigned char> color(common::taxonName2ID.nTaxa(), 0);
		for (int u = 0; u < n; ++u) {
			if (T.nd[u].taxon >= 0 || (int)T.nd[u].ch.size() != 2) continue;
			int const c0 = T.nd[u].ch[0], c1 = T.nd[u].ch[1];
			for (auto t : leaf[c0]) color[(size_t)t] = 1;   // M1 = 0b010
			for (auto t : leaf[c1]) color[(size_t)t] = 2;   // M2 = 0b100
			total += rscore::W_colors(g, color, sc);
			for (auto t : leaf[c0]) color[(size_t)t] = 0;
			for (auto t : leaf[c1]) color[(size_t)t] = 0;
		}
	}
	// score = Σ W / 2
	cout << (total / 2) << endl;
	return 0;
}
