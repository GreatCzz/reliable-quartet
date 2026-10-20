// solver: full-scale local search (NNI) for the reliable-quartet objective.
//
// score(S, genes) = sum over internal nodes v of S of F(child0(v), child1(v))
// where F(l, r) counts quartets ab|cd (a,b in l; c,d in r) that are reliable.
//
// An NNI at edge (u -> v) only changes the child sets of u and v, so only the
// two cached F values need recomputation -> O(2 * T * |g|) per candidate.
//
// Build:  g++ -std=c++20 -O3 -march=native -fopenmp -I ../aster2/src solver.cpp -o solver

#include "newick.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <thread>
#include <vector>
#include <omp.h>

namespace {

using std::size_t;
using std::string;
using std::vector;
using u64 = std::uint64_t;
using i64 = long long;

// ---------------------------------------------------------------- gene trees
struct GeneTree {
	vector<int> start;    // CSR offsets, size n+1
	vector<int> kids;     // CSR children
	vector<char> solid;   // edge (parent -> node) solid
	vector<int> taxon;    // leaf taxon, -1 internal
	int root = 0, n = 0;
};

GeneTree buildGene(newick::Newick::Node const* const root) {
	GeneTree g;
	vector<newick::Newick::Node const*> order, st{ root };
	while (!st.empty()) {
		auto* u = st.back(); st.pop_back();
		order.push_back(u);
		for (auto it = u->children.rbegin(); it != u->children.rend(); ++it) st.push_back(*it);
	}
	std::unordered_map<newick::Newick::Node const*, int> id;
	id.reserve(order.size() * 2);
	for (size_t i = 0; i < order.size(); ++i) id[order[i]] = (int)i;
	g.n = (int)order.size();
	g.start.assign(g.n + 1, 0);
	g.solid.assign(g.n, 1);
	g.taxon.assign(g.n, -1);
	for (size_t i = 0; i < order.size(); ++i)
		g.start[i + 1] = g.start[i] + (int)order[i]->children.size();
	g.kids.reserve(g.start[g.n]);
	for (size_t i = 0; i < order.size(); ++i) {
		auto* u = order[i];
		if (u->children.empty()) {
			g.taxon[i] = (int)common::taxonName2ID[u->name];
		}
		else {
			for (auto* c : u->children) {
				int const cid = id[c];
				g.kids.push_back(cid);
				g.solid[cid] = (c->name == "D") ? 0 : 1;
			}
		}
	}
	g.root = 0;
	return g;
}

// --------------------------------------------------------------- F evaluation
struct Scratch {
	vector<i64> c1, a1, m1, G1, H1, c2, a2, m2, G2, H2;
	void resize(int n) {
		c1.resize(n); a1.resize(n); m1.resize(n); G1.resize(n); H1.resize(n);
		c2.resize(n); a2.resize(n); m2.resize(n); G2.resize(n); H2.resize(n);
	}
};

inline u64 evalF(GeneTree const& g, u64 const* const bl, u64 const* const br, Scratch& s) {
	int const n = g.n;
	i64* const c1 = s.c1.data(); i64* const a1 = s.a1.data(); i64* const m1 = s.m1.data();
	i64* const c2 = s.c2.data(); i64* const a2 = s.a2.data(); i64* const m2 = s.m2.data();
	for (int i = n - 1; i >= 0; --i) {
		int const t = g.taxon[i];
		if (t >= 0) {
			c1[i] = ((bl[(size_t)t >> 6] >> (t & 63)) & 1ULL) ? 1 : 0;
			c2[i] = ((br[(size_t)t >> 6] >> (t & 63)) & 1ULL) ? 1 : 0;
			a1[i] = 0; m1[i] = 0; a2[i] = 0; m2[i] = 0;
		}
		else {
			i64 C1 = 0, X1 = 0, M1 = 0, C2 = 0, X2 = 0, M2 = 0;
			int const e = g.start[i + 1];
			for (int k = g.start[i]; k < e; ++k) {
				int const ch = g.kids[k];
				M1 += m1[ch]; M2 += m2[ch];
				if (!g.solid[ch]) continue;
				i64 const v1 = c1[ch], v2 = c2[ch];
				X1 += C1 * v1; C1 += v1;
				X2 += C2 * v2; C2 += v2;
			}
			c1[i] = C1; a1[i] = X1; m1[i] = X1 + M1;
			c2[i] = C2; a2[i] = X2; m2[i] = X2 + M2;
		}
	}
	i64 const sp = m1[g.root], sq = m2[g.root];
	if (sp == 0 || sq == 0) return 0;
	i64* const G1 = s.G1.data(); i64* const H1 = s.H1.data();
	i64* const G2 = s.G2.data(); i64* const H2 = s.H2.data();
	G1[g.root] = 0; H1[g.root] = 0; G2[g.root] = 0; H2[g.root] = 0;
	u64 corr = (u64)a1[g.root] * (u64)a2[g.root];
	for (int i = 0; i < n; ++i) {
		int const e = g.start[i + 1];
		for (int k = g.start[i]; k < e; ++k) {
			int const ch = g.kids[k];
			if (g.solid[ch]) {
				G1[ch] = (c1[i] - c1[ch]) + G1[i];
				G2[ch] = (c2[i] - c2[ch]) + G2[i];
			}
			else { G1[ch] = 0; G2[ch] = 0; }
			H1[ch] = c1[ch] * G1[ch];
			H2[ch] = c2[ch] * G2[ch];
			corr += (u64)a1[ch] * (u64)a2[ch] + (u64)a2[ch] * (u64)H1[ch] + (u64)a1[ch] * (u64)H2[ch];
		}
	}
	u64 const total = (u64)sp * (u64)sq;
	return total > corr ? total - corr : 0;
}

// --------------------------------------------------------------- species tree
struct SpeciesTree {
	int n = 0, root = 0;
	vector<int> par, c0, c1, taxon;
	vector<int> internal;      // node ids that have two children
	vector<u64> bits;          // per node, W words
	size_t W = 0;
};

SpeciesTree buildSpecies(newick::Newick::Node const* const root, size_t nTaxa) {
	SpeciesTree t;
	vector<newick::Newick::Node const*> order, st{ root };
	while (!st.empty()) {
		auto* u = st.back(); st.pop_back();
		order.push_back(u);
		for (auto it = u->children.rbegin(); it != u->children.rend(); ++it) st.push_back(*it);
	}
	std::unordered_map<newick::Newick::Node const*, int> id;
	id.reserve(order.size() * 2);
	for (size_t i = 0; i < order.size(); ++i) id[order[i]] = (int)i;
	t.n = (int)order.size();
	t.par.assign(t.n, -1); t.c0.assign(t.n, -1); t.c1.assign(t.n, -1); t.taxon.assign(t.n, -1);
	t.root = 0;
	t.W = (nTaxa + 63) / 64;
	t.bits.assign((size_t)t.n * t.W, 0);
	for (size_t i = 0; i < order.size(); ++i) {
		auto* u = order[i];
		if (u->children.empty()) {
			t.taxon[i] = (int)common::taxonName2ID[u->name];
		}
		else {
			if (u->children.size() != 2) {
				std::cerr << "Error: species tree must be binary" << std::endl;
				std::exit(1);
			}
			t.c0[i] = id[u->children[0]];
			t.c1[i] = id[u->children[1]];
			t.par[t.c0[i]] = (int)i;
			t.par[t.c1[i]] = (int)i;
			t.internal.push_back((int)i);
		}
	}
	for (int i = t.n - 1; i >= 0; --i) {
		u64* const dst = t.bits.data() + (size_t)i * t.W;
		if (t.taxon[i] >= 0) {
			size_t const x = (size_t)t.taxon[i];
			dst[x >> 6] |= 1ULL << (x & 63);
		}
		else {
			u64 const* const a = t.bits.data() + (size_t)t.c0[i] * t.W;
			u64 const* const b = t.bits.data() + (size_t)t.c1[i] * t.W;
			for (size_t w = 0; w < t.W; ++w) dst[w] = a[w] | b[w];
		}
	}
	return t;
}

// ----------------------------------------------------------------------- main
struct Args {
	string species, input, output;
	int threads = 8, sweeps = 20;
	int rounds = 1, kick = 0;
	bool verbose = false;
};

Args parseArgs(int argc, char** argv) {
	Args a;
	for (int i = 1; i < argc; ++i) {
		string s = argv[i];
		auto next = [&]() -> string { if (i + 1 < argc) return argv[++i]; std::exit(2); };
		if (s == "-s") a.species = next();
		else if (s == "-i") a.input = next();
		else if (s == "-o") a.output = next();
		else if (s == "-t") a.threads = std::stoi(next());
		else if (s == "--sweeps") a.sweeps = std::stoi(next());
		else if (s == "--rounds") a.rounds = std::stoi(next());
		else if (s == "--kick") a.kick = std::stoi(next());
		else if (s == "-v") a.verbose = true;
		else { std::cerr << "unknown arg " << s << std::endl; std::exit(2); }
	}
	return a;
}

string toNewick(SpeciesTree const& t, int u) {
	if (t.taxon[u] >= 0) return common::taxonName2ID[(size_t)t.taxon[u]];
	return "(" + toNewick(t, t.c0[u]) + "," + toNewick(t, t.c1[u]) + ")";
}

} // namespace

int main(int argc, char** argv) {
	Args const args = parseArgs(argc, argv);

	// species tree
	vector<newick::Newick> spTrees;
	{
		std::ifstream fin(args.species);
		string line;
		while (std::getline(fin, line)) {
			if (line.find_first_not_of(" \t\r\n") != string::npos) {
				spTrees.emplace_back(line);
				break;
			}
		}
	}
	vector<GeneTree> genes;
	{
		std::ifstream fin(args.input);
		string line;
		while (std::getline(fin, line)) {
			if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
			newick::Newick nw(line);
			genes.push_back(buildGene(nw.root()));
		}
	}
	size_t const nTaxa = common::taxonName2ID.nTaxa();
	SpeciesTree S = buildSpecies(spTrees[0].root(), nTaxa);

	int const T = (int)genes.size();
	int const K = (int)S.internal.size();
	std::cerr << "[solver] species nodes=" << S.n << " internal=" << K
		<< " genes=" << T << " taxa=" << nTaxa << std::endl;

	// cache: F(internal k, gene g)
	vector<i64> cache((size_t)K * (size_t)T, 0);

	int const P = std::max(1, args.threads);
	vector<Scratch> scratch(P);
	int maxN = 1;
	for (auto const& g : genes) maxN = std::max(maxN, g.n);
	for (auto& s : scratch) s.resize(maxN);

	auto t0 = std::chrono::steady_clock::now();
	#pragma omp parallel for schedule(dynamic) num_threads(P)
	for (int k = 0; k < K; ++k) {
		Scratch& sc = scratch[omp_get_thread_num()];
		int const v = S.internal[k];
		u64 const* const bl = S.bits.data() + (size_t)S.c0[v] * S.W;
		u64 const* const br = S.bits.data() + (size_t)S.c1[v] * S.W;
		for (int g = 0; g < T; ++g)
			cache[(size_t)k * T + g] = (i64)evalF(genes[g], bl, br, sc);
	}
	i64 total = 0;
	for (i64 x : cache) total += x;
	std::cerr << "[solver] initial score = " << total << "  ("
		<< std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - t0).count()
		<< " s initial scoring)" << std::endl;

	// node -> internal index
	vector<int> idx(S.n, -1);
	for (int k = 0; k < K; ++k) idx[S.internal[k]] = k;

	vector<u64> tmpA(S.W), tmpB(S.W);   // union scratch (per sweep, single threaded)

	auto doSweeps = [&](int maxSweeps) {
	for (int sweep = 1; sweep <= maxSweeps; ++sweep) {
		i64 sweepGain = 0;
		int moves = 0;
		auto ts = std::chrono::steady_clock::now();
		for (int k = 0; k < K; ++k) {
			int const v = S.internal[k];
			int const u = S.par[v];
			if (u < 0 || S.c0[u] < 0) continue;               // parent must be internal
			int const a = (S.c0[u] == v) ? S.c1[u] : S.c0[u];
			int const c = S.c0[v], d = S.c1[v];
			int const ku = idx[u];
			if (ku < 0) continue;

			// two NNI variants
			for (int var = 0; var < 2; ++var) {
				int const x = (var == 0) ? c : d;          // moves up to u
				int const y = (var == 0) ? d : c;          // stays with v
				// v's new children: (a, y); u's new children: (v, x)
				u64* const bA = S.bits.data() + (size_t)a * S.W;
				u64* const bY = S.bits.data() + (size_t)y * S.W;
				u64* const bX = S.bits.data() + (size_t)x * S.W;
				u64* const bV = S.bits.data() + (size_t)v * S.W;
				u64* const bU = S.bits.data() + (size_t)u * S.W;
				for (size_t w = 0; w < S.W; ++w) tmpA[w] = bA[w] | bY[w];   // new v
				// F(u_new) = F(tmpA, bX); F(v_new) = F(bA, bY)
				i64 delta = 0;
				#pragma omp parallel for schedule(dynamic) num_threads(P) reduction(+:delta)
				for (int g = 0; g < T; ++g) {
					Scratch& sc = scratch[omp_get_thread_num()];
					i64 const nu = (i64)evalF(genes[g], tmpA.data(), bX, sc);
					i64 const nv = (i64)evalF(genes[g], bA, bY, sc);
					delta += nu + nv - cache[(size_t)k * T + g] - cache[(size_t)ku * T + g];
				}
				if (delta > 0) {
					// accept: update tree + bits + cache
					S.c0[v] = a; S.c1[v] = y; S.par[a] = v; S.par[y] = v;
					S.c0[u] = v; S.c1[u] = x;
					S.par[v] = u; S.par[x] = u;
					memcpy(bV, tmpA.data(), S.W * sizeof(u64));
					// refresh cache for u and v
					for (int g = 0; g < T; ++g) {
						Scratch& sc = scratch[0];
						cache[(size_t)k * T + g] = (i64)evalF(genes[g], bA, bY, sc);
						cache[(size_t)ku * T + g] = (i64)evalF(genes[g], bV, bX, sc);
					}
					total += delta; sweepGain += delta; ++moves;
					break;   // move on to next edge
				}
			}
		}
		std::cerr << "[solver] sweep " << sweep << ": moves=" << moves << " gain=" << sweepGain
			<< " total=" << total << "  ("
			<< std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - ts).count()
			<< " s)" << std::endl;
		if (moves == 0) break;
	}
	};

	// ---- iterated local search: NNI to convergence, kick, repeat ----
	vector<int> bestC0 = S.c0, bestC1 = S.c1, bestPar = S.par;
	i64 bestTotal = total;
	for (int round = 1; round <= args.rounds; ++round) {
		if (round > 1) {
			// random kick
			std::mt19937_64 rng(12345 + round);
			for (int k = 0; k < args.kick; ++k) {
				// pick a random internal node v whose parent is internal
				int v = -1;
				for (int tries = 0; tries < 50; ++tries) {
					int const cand = S.internal[rng() % S.internal.size()];
					int const p = S.par[cand];
					if (p >= 0 && S.c0[p] >= 0) { v = cand; break; }
				}
				if (v < 0) break;
				int const u = S.par[v];
				int const a = (S.c0[u] == v) ? S.c1[u] : S.c0[u];
				int const c = S.c0[v], d = S.c1[v];
				int const ku = idx[u];
				int const x = (rng() & 1) ? c : d;
				int const y = (x == c) ? d : c;
				u64* const bA = S.bits.data() + (size_t)a * S.W;
				u64* const bY = S.bits.data() + (size_t)y * S.W;
				u64* const bX = S.bits.data() + (size_t)x * S.W;
				u64* const bV = S.bits.data() + (size_t)v * S.W;
				for (size_t w = 0; w < S.W; ++w) tmpA[w] = bA[w] | bY[w];
				S.c0[v] = a; S.c1[v] = y; S.par[a] = v; S.par[y] = v;
				S.c0[u] = v; S.c1[u] = x; S.par[v] = u; S.par[x] = u;
				memcpy(bV, tmpA.data(), S.W * sizeof(u64));
				for (int g = 0; g < T; ++g) {
					Scratch& sc = scratch[0];
					cache[(size_t)idx[v] * T + g] = (i64)evalF(genes[g], bA, bY, sc);
					cache[(size_t)ku * T + g] = (i64)evalF(genes[g], bV, bX, sc);
				}
				total = 0;
				for (i64 z : cache) total += z;
			}
		}
		doSweeps(args.sweeps);
		std::cerr << "[solver] round " << round << " total=" << total << std::endl;
		if (total > bestTotal) {
			bestTotal = total;
			bestC0 = S.c0; bestC1 = S.c1; bestPar = S.par;
		}
	}
	S.c0 = bestC0; S.c1 = bestC1; S.par = bestPar; total = bestTotal;

	string const nw = toNewick(S, S.root) + ";";
	if (args.output.empty()) std::cout << nw << std::endl;
	else { std::ofstream fout(args.output); fout << nw << "\n"; }
	std::cerr << "[solver] final score = " << total << std::endl;
	return 0;
}
