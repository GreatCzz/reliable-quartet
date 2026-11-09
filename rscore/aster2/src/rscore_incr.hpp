#ifndef RSCORE_INCR_HPP
#define RSCORE_INCR_HPP

// rscore_incr: 增量版可靠四叶树计分（13 状态 DP），接入 ASTER2 colorable 接口。
//
// 思路（见 rscore/FEASIBILITY_incremental.md）：
//  - 每个基因树先 Huffman resolve 成二叉（合成边标 C）；
//  - 维护两套 13 状态 DP（V1: C 当实线；V2: 只统计用到 C 的部分）；
//  - 改一个叶子的颜色只重算祖先链（O(深度)）；
//  - elementScore = V1[root][12] − V2[root][12]（= 2 × 节点分数，与 rscore.hpp 一致）。
//
// 注意：状态按"元素（基因）"存储；ASTER2 的 Task 按元素分片，同一元素同一时刻
// 只被一个线程处理，因此不需要每线程一份状态。

#include "stepwise_colorable.hpp"
#include "newick.hpp"

#include <array>
#include <algorithm>
#include <functional>
#include <queue>
#include <cstdint>
#include <cstdlib>
#include <unordered_map>
#include <vector>

namespace rscore_incr {

	using std::size_t;
	using std::string;
	using std::vector;
	using ll = long long;
	using St = std::array<ll, 13>;

	struct Attributes {
		using index_t = int;
		using score_t = long long;
		using support_t = bool;
		using length_t = bool;
		static inline bool constexpr WEIGHTING_BY_SUPPORT = false;
		static inline bool constexpr WEIGHTING_BY_LENGTH = false;
		static inline score_t constexpr ZERO = 0;
		static inline score_t constexpr EPSILON = 0;
	};

	struct GeneTree {
		vector<vector<int> > ch;
		vector<char> solid;
		vector<int> taxon;
		int root = -1;
		int n = 0;
	};

	inline GeneTree buildGeneTree(newick::Newick::Node const* const root) {
		GeneTree g;
		vector<newick::Newick::Node const*> order, st{ root };
		while (!st.empty()) {
			auto* u = st.back(); st.pop_back(); order.push_back(u);
			for (auto it = u->children.rbegin(); it != u->children.rend(); ++it) st.push_back(*it);
		}
		std::unordered_map<newick::Newick::Node const*, int> id;
		id.reserve(order.size() * 2);
		for (size_t i = 0; i < order.size(); ++i) id[order[i]] = (int)i;
		g.n = (int)order.size();
		g.ch.assign(g.n, {});
		g.solid.assign(g.n, 1);
		g.taxon.assign(g.n, -1);
		g.root = 0;
		for (size_t i = 0; i < order.size(); ++i) {
			auto* u = order[i];
			if (u->children.empty()) g.taxon[i] = (int)common::taxonName2ID[u->name];
			else for (auto* c : u->children) {
				int const cid = id[c];
				g.ch[i].push_back(cid);
				g.solid[cid] = (c->name == "D") ? 0 : 1;
			}
		}
		return g;
	}

	// ---------------- 二叉化 + 13 状态 ----------------
	// code: 0=C(合成) 1=S 2=D 3=leaf
	struct BinaryGene {
		vector<int> lc, rc, par, order, sz;
		vector<unsigned char> code;
		vector<int> leafTaxon;           // 每节点：叶的 taxon id，内部为 -1
		vector<vector<St> > V;           // [2][node] -> 13 状态（V[0]=V1, V[1]=V2）
		std::unordered_map<int, int> taxonNode;
		int root = 0, n = 0;
	};

	inline void step(St& A, const St& U, const St& V, int cu, int cv, int nu, int nv, bool second) {
		ll Ru=U[0],Bu=U[1],Gu=U[2],RRu=U[3],BBu=U[4],GGu=U[5],RBu=U[6],RGu=U[7],BGu=U[8],coRu=U[9],coBu=U[10],coGu=U[11],Qu=U[12];
		ll Rv=V[0],Bv=V[1],Gv=V[2],RRv=V[3],BBv=V[4],GGv=V[5],RBv=V[6],RGv=V[7],BGv=V[8],coRv=V[9],coBv=V[10],coGv=V[11],Qv=V[12];
		A[0]=nu*Ru+nv*Rv; A[1]=nu*Bu+nv*Bv; A[2]=nu*Gu+nv*Gv;
		if (!second) {
			A[3]=RRu+RRv+nu*nv*Ru*Rv; A[4]=BBu+BBv+nu*nv*Bu*Bv; A[5]=GGu+GGv+nu*nv*Gu*Gv;
			A[6]=RBu+RBv+nu*nv*(Ru*Bv+Bu*Rv);
			A[7]=RGu+RGv+nu*nv*(Ru*Gv+Gu*Rv);
			A[8]=BGu+BGv+nu*nv*(Bu*Gv+Gu*Bv);
			A[9]=nu*coRu+nv*coRv + (BBu*nv*Gv + BBv*nu*Gu + GGu*nv*Bv + GGv*nu*Bu + BGu*nv*Rv + BGv*nu*Ru);
			A[10]=nu*coBu+nv*coBv + (RRu*nv*Gv + RRv*nu*Gu + GGu*nv*Rv + GGv*nu*Ru + RGu*nv*Bv + RGv*nu*Bu);
			A[11]=nu*coGu+nv*coGv + (BBu*nv*Rv + BBv*nu*Ru + RRu*nv*Bv + RRv*nu*Bu + RBu*nv*Gv + RBv*nu*Gu);
			A[12]=Qu+Qv + nu*nv*(coRu*Rv+coRv*Ru+coBu*Bv+coBv*Bu+coGu*Gv+coGv*Gu)
				+ (RRu*BGv+RRv*BGu+BBu*RGv+BBv*RGu+GGu*RBv+GGv*RBu);
		} else {
			A[3]=cu*RRu+cv*RRv+nu*nv*Ru*Rv; A[4]=cu*BBu+cv*BBv+nu*nv*Bu*Bv; A[5]=cu*GGu+cv*GGv+nu*nv*Gu*Gv;
			A[6]=cu*RBu+cv*RBv+nu*nv*(Ru*Bv+Bu*Rv);
			A[7]=cu*RGu+cv*RGv+nu*nv*(Ru*Gv+Gu*Rv);
			A[8]=cu*BGu+cv*BGv+nu*nv*(Bu*Gv+Gu*Bv);
			A[9]=nu*coRu+nv*coRv + cu*nv*(BBu*Gv+GGu*Bv+BGu*Rv) + cv*nu*(BBv*Gu+GGv*Bu+BGv*Ru);
			A[10]=nu*coBu+nv*coBv + cu*nv*(RRu*Gv+GGu*Rv+RGu*Bv) + cv*nu*(RRv*Gu+GGv*Ru+RGv*Bu);
			A[11]=nu*coGu+nv*coGv + cu*nv*(RRu*Bv+BBu*Rv+RBu*Gv) + cv*nu*(RRv*Bu+BBv*Ru+RBv*Gu);
			A[12]=Qu+Qv + nu*nv*(coRu*Rv+coRv*Ru+coBu*Bv+coBv*Bu+coGu*Gv+coGv*Gu)
				+ cu*cv*(RRu*BGv+RRv*BGu+BBu*RGv+BBv*RGu+GGu*RBv+GGv*RBu);
		}
	}

	inline BinaryGene makeBinary(GeneTree const& g) {
		BinaryGene b;
		std::function<int(int)> conv = [&](int u) -> int {
			int id = b.n++;
			b.lc.push_back(-1); b.rc.push_back(-1); b.par.push_back(-1); b.sz.push_back(0);
			b.code.push_back(g.ch[u].empty() ? 3 : (g.solid[u] ? 1 : 2));
			b.leafTaxon.push_back(g.taxon[u]);
			if (g.ch[u].empty()) { b.sz[id] = 1; return id; }
			vector<int> ids;
			for (int c : g.ch[u]) ids.push_back(conv(c));
			if (ids.size() == 1) {
				int c = ids[0];
				b.lc[id] = b.lc[c]; b.rc[id] = b.rc[c]; b.code[id] = b.code[c];
				b.sz[id] = b.sz[c]; b.leafTaxon[id] = b.leafTaxon[c];
				if (b.lc[id] >= 0) { b.par[b.lc[id]] = id; b.par[b.rc[id]] = id; }
				return id;
			}
			if (ids.size() == 2) { b.lc[id] = ids[0]; b.rc[id] = ids[1]; }
			else {
				std::priority_queue<std::pair<int,int>, vector<std::pair<int,int>>, std::greater<>> pq;
				for (int c : ids) pq.push({ b.sz[c], c });
				while (pq.size() > 2) {
					auto [s1, a] = pq.top(); pq.pop();
					auto [s2, c] = pq.top(); pq.pop();
					int nid = b.n++;
					b.lc.push_back(a); b.rc.push_back(c); b.par.push_back(-1);
					b.sz.push_back(s1 + s2); b.code.push_back(0); b.leafTaxon.push_back(-1);
					b.par[a] = nid; b.par[c] = nid;
					pq.push({ s1 + s2, nid });
				}
				auto [s1, a] = pq.top(); pq.pop();
				auto [s2, c] = pq.top(); pq.pop();
				b.lc[id] = a; b.rc[id] = c;
			}
			b.par[b.lc[id]] = id; b.par[b.rc[id]] = id;
			b.sz[id] = b.sz[b.lc[id]] + b.sz[b.rc[id]];
			return id;
		};
		b.root = conv(g.root);
		vector<int> st{ b.root };
		while (!st.empty()) {
			int t = st.back(); st.pop_back();
			b.order.push_back(t);
			if (b.lc[t] >= 0) { st.push_back(b.lc[t]); st.push_back(b.rc[t]); }
		}
		std::reverse(b.order.begin(), b.order.end());
		for (int i = 0; i < b.n; ++i)
			if (b.lc[i] < 0 && b.leafTaxon[i] >= 0) b.taxonNode[b.leafTaxon[i]] = i;
		return b;
	}

	inline void initStates(BinaryGene& b) {
		b.V.assign(2, vector<St>(b.n));
		for (auto& v : b.V) for (auto& s : v) s.fill(0);
		for (int t : b.order) {
			if (b.lc[t] < 0) { b.V[0][t][0] = 1; b.V[1][t][0] = 1; continue; }
			int l = b.lc[t], r = b.rc[t];
			int nu = (b.code[l] == 2) ? 0 : 1, nv = (b.code[r] == 2) ? 0 : 1;
			int cu = (b.code[l] == 0) ? 1 : 0, cv = (b.code[r] == 0) ? 1 : 0;
			step(b.V[0][t], b.V[0][l], b.V[0][r], cu, cv, nu, nv, false);
			step(b.V[1][t], b.V[1][l], b.V[1][r], cu, cv, nu, nv, true);
		}
	}

	inline void recolor(BinaryGene& b, int leaf, int c) {
		b.V[0][leaf].fill(0); b.V[1][leaf].fill(0);
		b.V[0][leaf][c] = 1; b.V[1][leaf][c] = 1;
		for (int a = b.par[leaf]; a >= 0; a = b.par[a]) {
			int l = b.lc[a], r = b.rc[a];
			int nu = (b.code[l] == 2) ? 0 : 1, nv = (b.code[r] == 2) ? 0 : 1;
			int cu = (b.code[l] == 0) ? 1 : 0, cv = (b.code[r] == 0) ? 1 : 0;
			step(b.V[0][a], b.V[0][l], b.V[0][r], cu, cv, nu, nv, false);
			step(b.V[1][a], b.V[1][l], b.V[1][r], cu, cv, nu, nv, true);
		}
	}

	// ---------------- 四色 2-2 计数（NNI 用） ----------------
	// 直接复用可靠的 F(X,Y)（两趟 DP，支持多叉，无需 resolve）。
	// 颜色 0..3；三种配对：{0,1}|{2,3}、{0,2}|{1,3}、{0,3}|{1,2}。
	struct FScratch {
		vector<ll> c1, a1, m1, G1, H1, c2, a2, m2, G2, H2;
		void resize(int n) {
			c1.resize(n); a1.resize(n); m1.resize(n); G1.resize(n); H1.resize(n);
			c2.resize(n); a2.resize(n); m2.resize(n); G2.resize(n); H2.resize(n);
		}
	};

	inline unsigned long long F_colors(GeneTree const& g, vector<unsigned char> const& color,
		int maskX, int maskY, FScratch& s) {
		int const n = g.n;
		ll* c1 = s.c1.data(); ll* a1 = s.a1.data(); ll* m1 = s.m1.data();
		ll* c2 = s.c2.data(); ll* a2 = s.a2.data(); ll* m2 = s.m2.data();
		for (int i = n - 1; i >= 0; --i) {
			int const t = g.taxon[i];
			if (t >= 0) {
				unsigned char const col = color[(size_t)t];
				c1[i] = ((maskX >> col) & 1) ? 1 : 0;
				c2[i] = ((maskY >> col) & 1) ? 1 : 0;
				a1[i] = 0; m1[i] = 0; a2[i] = 0; m2[i] = 0;
			} else {
				ll C1 = 0, X1 = 0, M1 = 0, C2 = 0, X2 = 0, M2 = 0;
				for (int c : g.ch[i]) {
					M1 += m1[c]; M2 += m2[c];
					if (!g.solid[c]) continue;
					ll const v1 = c1[c], v2 = c2[c];
					X1 += C1 * v1; C1 += v1;
					X2 += C2 * v2; C2 += v2;
				}
				c1[i] = C1; a1[i] = X1; m1[i] = X1 + M1;
				c2[i] = C2; a2[i] = X2; m2[i] = X2 + M2;
			}
		}
		ll const sp = m1[g.root], sq = m2[g.root];
		if (sp == 0 || sq == 0) return 0;
		ll* G1 = s.G1.data(); ll* H1 = s.H1.data(); ll* G2 = s.G2.data(); ll* H2 = s.H2.data();
		G1[g.root] = 0; H1[g.root] = 0; G2[g.root] = 0; H2[g.root] = 0;
		unsigned long long corr = (unsigned long long)a1[g.root] * (unsigned long long)a2[g.root];
		for (int i = 0; i < n; ++i) for (int c : g.ch[i]) {
			if (g.solid[c]) { G1[c] = (c1[i]-c1[c]) + G1[i]; G2[c] = (c2[i]-c2[c]) + G2[i]; }
			else { G1[c] = 0; G2[c] = 0; }
			H1[c] = c1[c] * G1[c]; H2[c] = c2[c] * G2[c];
			corr += (unsigned long long)a1[c]*(unsigned long long)a2[c]
				+ (unsigned long long)a2[c]*(unsigned long long)H1[c]
				+ (unsigned long long)a1[c]*(unsigned long long)H2[c];
		}
		unsigned long long const total = (unsigned long long)sp * (unsigned long long)sq;
		return total > corr ? total - corr : 0;
	}

	// ---------------- ASTER2 colorable 接口 ----------------
	template<class Attributes> class Color {
	public:
		using score_t = typename Attributes::score_t;
		using index_t = typename Attributes::index_t;
		static inline bool constexpr IS_ROOTED = false;
		static inline score_t constexpr ZERO = Attributes::ZERO;
		static inline score_t constexpr EPSILON = Attributes::EPSILON;

		struct SharedConstData {
			using ParentClass = Color<Attributes>;
			vector<GeneTree> genes;
			size_t nTaxa = 0;
			size_t nElements() const noexcept { return genes.size(); }
		};

		explicit Color(SharedConstData const& data) noexcept : data_(&data), trees_(data.genes.size(), BinaryGene()),
			quadColor_(data.genes.size(), vector<unsigned char>(data.nTaxa, 0)) {
			size_t const maxGenes = [&]() -> size_t {
				const char* e = std::getenv("RSCORE_MAXGENES");
				return e ? (size_t)std::atoi(e) : data.genes.size();
			}();
			size_t const k = std::min(maxGenes, data.genes.size());
			for (size_t i = 0; i < k; ++i) {
				trees_[i] = makeBinary(data.genes[i]);
				initStates(trees_[i]);
			}
		}

		void elementSetTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept {
			if (iElement < quadColor_.size()) {
				if (iTaxon < data_->nTaxa) quadColor_[iElement][iTaxon] = (unsigned char)iColor;
				if (quadMode_) return;   // 四色模式：只更新颜色，3 色 DP 不动
			}
			if (iElement >= trees_.size()) return;
			auto& b = trees_[iElement];
			if (b.n == 0) return;
			auto it = b.taxonNode.find((int)iTaxon);
			if (it == b.taxonNode.end()) return;   // 该基因树缺这个类群
			recolor(b, it->second, (int)iColor);
		}

		void elementClearTaxonColor(size_t iElement, size_t iTaxon, size_t) noexcept {
			elementSetTaxonColor(iElement, iTaxon, 0);
		}

		void elementClearAndSetTaxonColor(size_t iElement, size_t iTaxon, size_t, size_t jColor) noexcept {
			elementSetTaxonColor(iElement, iTaxon, jColor);
		}

		score_t elementScore(size_t iElement) noexcept {
			if (iElement >= trees_.size()) return ZERO;
			auto const& b = trees_[iElement];
			if (b.n == 0) return ZERO;
			return (score_t)(b.V[0][b.root][12] - b.V[1][b.root][12]);
		}

	private:
		SharedConstData const* data_;
		vector<BinaryGene> trees_;
		vector<vector<unsigned char> > quadColor_;   // 四色模式下每元素每类群的颜色
		bool quadMode_ = false;
	};

	class Documentation : public common::DocumentationBase {
	protected:
		string introduction() const noexcept override {
			return R"DOC(# rscore_incr: incremental reliable-quartet species tree
)DOC";
		}
		string input() const noexcept override {
			return R"DOC(# INPUT
One gene tree per line, Newick format. Internal edge labels "S"/"D" mark solid
and dashed edges; leaf edges are always solid. Polytomies are allowed.
)DOC";
		}
		string exampleInput() const noexcept override {
			return R"DOC(# EXAMPLE INPUT
```
((a,b)S,(c,d)D);
((a,b)S,c,d);
```
)DOC";
		}
	public:
		string programName() const noexcept override { return "rscore_incr"; }
	};

	template<bool> class Driver : public common::LogInfo {
	public:
		using DataClasses = std::variant<typename Color<Attributes>::SharedConstData>;

		static std::pair<string, string> programNames() noexcept {
			return { "rscore_incr", "Incremental reliable quartet score" };
		}

		static void addArguments() noexcept {}

		static DataClasses getStepwiseColorSharedConstData() noexcept {
			return read();
		}

	private:
		static typename Color<Attributes>::SharedConstData read() noexcept {
			typename Color<Attributes>::SharedConstData data;
			std::ifstream fin(ARG.get<string>("input"));
			if (!fin.is_open()) {
				std::cerr << "Error: unable to open input file" << std::endl;
				std::exit(-1);
			}
			string line;
			while (std::getline(fin, line)) {
				if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
				newick::Newick nw(line);
				data.genes.push_back(buildGeneTree(nw.root()));
			}
			data.nTaxa = common::taxonName2ID.nTaxa();
			return data;
		}
	};

} // namespace rscore_incr

#endif
