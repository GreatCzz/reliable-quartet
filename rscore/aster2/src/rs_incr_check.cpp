// rs_incr_check: 用 rscore_incr 的增量 Color 算一个三分划的 W，和 rscore::W_colors 对拍，并测改色成本
#include "rscore_incr.hpp"
#include "rscore.hpp"
#include <chrono>
#include <fstream>
#include <iostream>
using namespace std;

int main(int argc, char** argv) {
	string spFile = argv[1], gFile = argv[2];
	int nG = argc > 3 ? atoi(argv[3]) : 20;
	newick::Newick sp;
	{ ifstream fin(spFile); string line; while (getline(fin, line)) if (line.find_first_not_of(" \t\r\n") != string::npos) { sp = newick::Newick(line); break; } }
	// 物种树：取根的左右孩子做三分划 A/B
	// 三分划：A = root->left->left 的叶子，B = root->left->right 的叶子，rest = root->right
	vector<int> A, B;
	{
		auto* L = sp.root()->children[0];
		auto* R = sp.root()->children[1];
		auto collect = [&](newick::Newick::Node const* u, vector<int>& out) {
			vector<newick::Newick::Node const*> st{u};
			while (!st.empty()) { auto* x = st.back(); st.pop_back();
				if (x->children.empty()) out.push_back((int)common::taxonName2ID[x->name]);
				else for (auto* c : x->children) st.push_back(c); }
		};
		vector<int> all;
		collect(L, all);
		if (L->children.size() >= 2) { collect(L->children[0], A); collect(L->children[1], B); }
		else { A = all; }
		(void)R;
	}
	rscore_incr::Color<rscore_incr::Attributes>::SharedConstData data;
	data.nTaxa = common::taxonName2ID.nTaxa();
	{ ifstream fin(gFile); string line; int cnt = 0;
	  while (getline(fin, line) && cnt < nG) {
		if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
		newick::Newick nw(line); data.genes.push_back(rscore_incr::buildGeneTree(nw.root())); ++cnt; } }
	rscore_incr::Color<rscore_incr::Attributes> col(data);

	// 设置三分划：A -> 1, B -> 2, 其余 0
	auto t0 = chrono::steady_clock::now();
	for (size_t e = 0; e < data.genes.size(); ++e) {
		for (int t : A) col.elementSetTaxonColor(e, (size_t)t, 1);
		for (int t : B) col.elementSetTaxonColor(e, (size_t)t, 2);
	}
	long long incr = 0;
	for (size_t e = 0; e < data.genes.size(); ++e) incr += col.elementScore(e);
	auto t1 = chrono::steady_clock::now();

	// 参照：rscore::W_colors（用另一套 GeneTree 类型，重新解析）
	vector<rscore::GeneTree> refGenes;
	{ ifstream fin(gFile); string line; int cnt = 0;
	  while (getline(fin, line) && cnt < nG) {
		if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
		newick::Newick nw(line); refGenes.push_back(rscore::buildGeneTree(nw.root())); ++cnt; } }
	long long ref = 0;
	for (auto const& g : refGenes) {
		rscore::Scratch sc; sc.resize(g.n);
		vector<unsigned char> color(common::taxonName2ID.nTaxa(), 0);
		for (int t : A) color[(size_t)t] = 1;
		for (int t : B) color[(size_t)t] = 2;
		ref += rscore::W_colors(g, color, sc);
	}
	double const sec = chrono::duration<double>(t1 - t0).count();
	cerr << "[incr] genes=" << data.genes.size() << " A=" << A.size() << " B=" << B.size()
	     << "\n[incr] elementScore sum = " << incr << "   W_colors sum = " << ref
	     << "   match=" << (incr == ref ? "YES" : "NO")
	     << "\n[incr] set+read time = " << sec << " s  (" << (sec / (double)data.genes.size() * 1e6) << " us/gene)\n";
	return incr == ref ? 0 : 1;
}
