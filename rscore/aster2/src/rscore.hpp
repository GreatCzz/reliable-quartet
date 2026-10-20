#ifndef RSCORE_HPP
#define RSCORE_HPP

// rscore: species-tree inference under the "reliable quartet" objective.
//
// This tool plugs the reliable-quartet score into ASTER2's stepwise-colouring
// iteration framework: the framework colours, for every internal node of the
// candidate species tree, the taxa of the two child clusters (colour 2 and 1)
// and the rest (colour 0).  Our objective decomposes exactly over those nodes:
//
//     score = sum over internal nodes v of F(cluster(colour 2), cluster(colour 1))
//
// where F(l, r) counts quartets ab|cd with a,b in l, c,d in r that are
// "reliable" in the gene tree (paths a->b and c->d are entirely solid).
// Hence elementScore() == F, and the framework's accumulated tree score is
// exactly the official reliable-quartet score.
//
// R1: non-incremental implementation (O(n) per elementScore call).

#include "stepwise_colorable.hpp"
#include "newick.hpp"

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <vector>

namespace rscore {

	using std::size_t;
	using std::string;
	using std::vector;

	struct Attributes {
		using index_t = int;
		using score_t = unsigned long long;
		using support_t = bool;
		using length_t = bool;
		static inline bool constexpr WEIGHTING_BY_SUPPORT = false;
		static inline bool constexpr WEIGHTING_BY_LENGTH = false;
		static inline score_t constexpr ZERO = 0;
		static inline score_t constexpr EPSILON = 0;
	};

	struct GeneTree {
		vector<vector<int> > ch;
		vector<char> solid;   // solid[c] = edge (parent(c) -> c) is solid
		vector<int> taxon;    // leaf taxon id, -1 for internal nodes
		int root = -1;
		int n = 0;
	};

	inline GeneTree buildGeneTree(newick::Newick::Node const* const root) {
		GeneTree g;
		vector<newick::Newick::Node const*> order;
		vector<newick::Newick::Node const*> st;
		st.push_back(root);
		while (!st.empty()) {
			newick::Newick::Node const* const u = st.back();
			st.pop_back();
			order.push_back(u);
			for (auto it = u->children.rbegin(); it != u->children.rend(); ++it)
				st.push_back(*it);
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
			newick::Newick::Node const* const u = order[i];
			if (u->children.empty()) {
				g.taxon[i] = (int)common::taxonName2ID[u->name];
			}
			else {
				for (newick::Newick::Node const* const c : u->children) {
					int const cid = id[c];
					g.ch[i].push_back(cid);
					g.solid[cid] = (c->name == "D") ? 0 : 1;
				}
			}
		}
		return g;
	}

	// Scratch buffers for one reliable-quartet evaluation.
	struct Scratch {
		vector<long long> c1, a1, m1, G1, H1;
		vector<long long> c2, a2, m2, G2, H2;
		void resize(int n) {
			c1.resize(n); a1.resize(n); m1.resize(n); G1.resize(n); H1.resize(n);
			c2.resize(n); a2.resize(n); m2.resize(n); G2.resize(n); H2.resize(n);
		}
	};

	// F(cluster of colour 1, cluster of colour 2) for one gene tree.
	inline std::uint64_t reliableScore(GeneTree const& g,
		vector<unsigned char> const& color, Scratch& s) {
		int const n = g.n;
		long long* const c1 = s.c1.data();
		long long* const a1 = s.a1.data();
		long long* const m1 = s.m1.data();
		long long* const c2 = s.c2.data();
		long long* const a2 = s.a2.data();
		long long* const m2 = s.m2.data();

		// bottom-up: cnt / anc (pairs anchored here) / cum (all pairs in subtree)
		for (int i = n - 1; i >= 0; --i) {
			int const t = g.taxon[i];
			if (t >= 0) {
				unsigned char const col = color[(size_t)t];
				c1[i] = (col == 1) ? 1 : 0; a1[i] = 0; m1[i] = 0;
				c2[i] = (col == 2) ? 1 : 0; a2[i] = 0; m2[i] = 0;
			}
			else {
				long long C1 = 0, X1 = 0, M1 = 0, C2 = 0, X2 = 0, M2 = 0;
				for (int const ch : g.ch[i]) {
					M1 += m1[ch];          // 子树内部的成对与父边无关，必须无条件累加
					M2 += m2[ch];
					if (!g.solid[ch]) continue;
					long long const v1 = c1[ch], v2 = c2[ch];
					X1 += C1 * v1; C1 += v1;
					X2 += C2 * v2; C2 += v2;
				}
				c1[i] = C1; a1[i] = X1; m1[i] = X1 + M1;
				c2[i] = C2; a2[i] = X2; m2[i] = X2 + M2;
			}
		}
		long long const sp = m1[g.root];
		long long const sq = m2[g.root];
		if (sp == 0 || sq == 0) return 0;

		// top-down: G / H and the three correction sums
		long long* const G1 = s.G1.data();
		long long* const H1 = s.H1.data();
		long long* const G2 = s.G2.data();
		long long* const H2 = s.H2.data();
		G1[g.root] = 0; H1[g.root] = 0;
		G2[g.root] = 0; H2[g.root] = 0;
		std::uint64_t corr = (std::uint64_t)a1[g.root] * (std::uint64_t)a2[g.root];
		for (int i = 0; i < n; ++i) {
			for (int const ch : g.ch[i]) {
				if (g.solid[ch]) {
					G1[ch] = (c1[i] - c1[ch]) + G1[i];
					G2[ch] = (c2[i] - c2[ch]) + G2[i];
				}
				else {
					G1[ch] = 0;
					G2[ch] = 0;
				}
				H1[ch] = c1[ch] * G1[ch];
				H2[ch] = c2[ch] * G2[ch];
				corr += (std::uint64_t)a1[ch] * (std::uint64_t)a2[ch]
					+ (std::uint64_t)a2[ch] * (std::uint64_t)H1[ch]
					+ (std::uint64_t)a1[ch] * (std::uint64_t)H2[ch];
			}
		}
		std::uint64_t const total = (std::uint64_t)sp * (std::uint64_t)sq;
		if (corr > total) {
			static int warned = 0;
			if (warned < 3) {
				++warned;
				std::cerr << "[rscore] corr > total: corr=" << corr << " total=" << total
					<< " sp=" << sp << " sq=" << sq << std::endl;
			}
			return 0;
		}
		return total - corr;
	}

	template<class Attributes> class Color {
	public:
		using score_t = typename Attributes::score_t;
		using support_t = typename Attributes::support_t;
		using length_t = typename Attributes::length_t;
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

		struct Tree {
			vector<unsigned char> color;   // per taxon: 0 rest, 1 sister, 2 subtree
			explicit Tree(size_t nTaxa) : color(nTaxa, 0) {}
		};

		explicit Color(SharedConstData const& data) noexcept
			: data_(&data), trees_(data.genes.size(), Tree(data.nTaxa)) {
			int maxN = 1;
			for (auto const& g : data.genes) if (g.n > maxN) maxN = g.n;
			scratch_.resize(maxN);
		}

		void elementSetTaxonColor(size_t iElement, size_t iTaxon, size_t iColor) noexcept {
			trees_[iElement].color[iTaxon] = (unsigned char)iColor;
		}

		void elementClearTaxonColor(size_t iElement, size_t iTaxon, size_t) noexcept {
			trees_[iElement].color[iTaxon] = 0;
		}

		void elementClearAndSetTaxonColor(size_t iElement, size_t iTaxon, size_t, size_t jColor) noexcept {
			trees_[iElement].color[iTaxon] = (unsigned char)jColor;
		}

		score_t elementScore(size_t iElement) noexcept {
			static thread_local Scratch tl;
			int maxN = 1;
			if ((int)tl.c1.size() < data_->genes[iElement].n) { tl.resize(data_->genes[iElement].n); }
			auto const& col = trees_[iElement].color;
			score_t const res = (score_t)reliableScore(data_->genes[iElement], col, tl);
			static std::atomic<int> nDbg{ 0 };
			if (std::getenv("RSCORE_DEBUG") && nDbg.fetch_add(1) < 10) {
				size_t n1 = 0, n2 = 0;
				for (auto c : col) { if (c == 1) ++n1; else if (c == 2) ++n2; }
				std::cerr << "[dbg] elem=" << iElement << " |color1|=" << n1
					<< " |color2|=" << n2 << " score=" << res << std::endl;
			}
			return res;
		}

	private:
		SharedConstData const* data_;
		vector<Tree> trees_;
		Scratch scratch_;   // NOTE: R1 is single-threaded per element only if called from one thread
	};

	class Documentation : public common::DocumentationBase {
	protected:
		string introduction() const noexcept override {
			return R"DOC(# rscore: reliable-quartet species tree (baseline)
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
		string programName() const noexcept override { return "rscore"; }
	};

	template<bool> class Driver : public common::LogInfo {
	public:
		using DataClasses = std::variant<typename Color<Attributes>::SharedConstData>;

		static std::pair<string, string> programNames() noexcept {
			return { "rscore", "Reliable quartet score (baseline)" };
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

} // namespace rscore

#endif // !RSCORE_HPP
