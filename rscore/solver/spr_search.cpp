// Two-stage restricted SPR search harness.
// Experimental only: it does not alter builder/solver defaults.
#define SPR_ORACLE_NO_MAIN
#include "spr_oracle.cpp"

#include <array>
#include <set>
#include <tuple>

namespace {

struct MoveKey {
    int p, r, q, t;
    bool operator<(MoveKey const& o) const {
        return std::tie(p, r, q, t) < std::tie(o.p, o.r, o.q, o.t);
    }
};

MoveKey keyOf(SPRMove const& m) { return {m.cutParent, m.cutChild, m.targetParent, m.targetChild}; }

int cutBin(int size) {
    if (size <= 1) return 0;
    if (size == 2) return 1;
    if (size <= 4) return 2;
    if (size <= 8) return 3;
    if (size <= 16) return 4;
    if (size <= 32) return 5;
    if (size <= 64) return 6;
    return 7;
}

char const* cutBinName(int b) {
    static char const* names[] = {"1", "2", "3-4", "5-8", "9-16", "17-32", "33-64", "65+"};
    return names[b];
}

vector<SPRMove> selectCuts(SpeciesTree const& t, int cutsPerBin, int maxCuts,
                          std::uint64_t seed) {
    std::array<vector<SPRMove>, 8> buckets;
    for (int p : t.internal) {
        if (t.par[p] < 0) continue;
        for (int r : {t.c0[p], t.c1[p]})
            buckets[(size_t)cutBin(t.sz[r])].push_back({p, r, -1, -1});
    }
    std::mt19937_64 rng(seed);
    vector<SPRMove> out;
    for (int b = 0; b < 8; ++b) {
        auto& v = buckets[(size_t)b];
        std::shuffle(v.begin(), v.end(), rng);
        int const take = std::min<int>(cutsPerBin, v.size());
        for (int i = 0; i < take; ++i) {
            out.push_back(v[(size_t)i]);
            std::cerr << "[spr-search] cut-bin " << cutBinName(b)
                      << " size=" << t.sz[v[(size_t)i].cutChild] << '\n';
        }
    }
    if (maxCuts > 0 && (int)out.size() > maxCuts) {
        std::shuffle(out.begin(), out.end(), rng);
        out.resize((size_t)maxCuts);
    }
    return out;
}

vector<SPRMove> allMovesForCut(SpeciesTree const& t, int p, int r, int maxTargets,
                               std::uint64_t seed) {
    vector<SPRMove> out;
    for (int q : t.internal) for (int target : {t.c0[q], t.c1[q]}) {
        SPRMove m{p, r, q, target};
        if (validMove(t, m)) out.push_back(m);
    }
    if (maxTargets > 0 && (int)out.size() > maxTargets) {
        std::mt19937_64 rng(seed);
        std::shuffle(out.begin(), out.end(), rng);
        out.resize((size_t)maxTargets);
    }
    return out;
}

struct ScoredMove {
    SPRMove move;
    int bucket = -1;
    int cutSize = 0;
    __int128 delta = 0;
};

bool validateCandidate(SpeciesTree const& base, SPRMove const& m,
                       vector<GeneTree> const& genes, int T, int threads,
                       vector<i64> const& oldCache, __int128 expected) {
    SpeciesTree clone = base;
    applySPR(clone, m);
    int const nTaxa = base.sz[base.root];
    if (!checkInvariants(clone, nTaxa) || clone.internal != base.internal) return false;
    auto fresh = refreshCache(clone, genes, T, threads);
    __int128 const actual = sumCache(fresh) - sumCache(oldCache);
    return actual == expected && (actual & 1) == 0;
}

} // namespace

int main(int argc, char** argv) {
    string speciesFile, genesFile, outputFile;
    int threads = 8, screenGenes = 8, cutsPerBin = 8, maxCuts = 0, maxTargets = 0;
    int perBin = 4, globalTop = 16;
    std::uint64_t seed = 20261111;
    for (int i = 1; i < argc; ++i) {
        string s = argv[i];
        auto next = [&]() -> string { if (++i >= argc) std::exit(2); return argv[i]; };
        if (s == "-s") speciesFile = next();
        else if (s == "-i") genesFile = next();
        else if (s == "-o") outputFile = next();
        else if (s == "-t") threads = std::stoi(next());
        else if (s == "--screen-genes") screenGenes = std::stoi(next());
        else if (s == "--cuts-per-bin") cutsPerBin = std::stoi(next());
        else if (s == "--max-cuts") maxCuts = std::stoi(next());
        else if (s == "--max-targets") maxTargets = std::stoi(next());
        else if (s == "--per-bin") perBin = std::stoi(next());
        else if (s == "--global-top") globalTop = std::stoi(next());
        else if (s == "--seed") seed = std::stoull(next());
        else { std::cerr << "unknown arg " << s << '\n'; return 2; }
    }
    if (speciesFile.empty() || genesFile.empty() || threads < 1 || screenGenes < 1
        || cutsPerBin < 1 || maxCuts < 0 || perBin < 0 || globalTop < 1) return 2;

    auto expectedNames = readTaxonNames(genesFile);
    auto genes = readGenes(genesFile);
    size_t const nTaxa = common::taxonName2ID.nTaxa();
    auto base = readSpecies(speciesFile, nTaxa, expectedNames);
    if (!checkInvariants(base, (int)nTaxa)) return 1;
    int const Tscreen = std::min<int>(screenGenes, genes.size());
    int const Tfull = (int)genes.size();
    auto oldScreen = refreshCache(base, genes, Tscreen, threads);
    auto oldFull = refreshCache(base, genes, Tfull, threads);
    std::cerr << "[spr-search] baseline score="
              << printInt128(sumCache(oldFull) / 2) << " taxa=" << nTaxa
              << " genes=" << Tfull << '\n';

    auto cuts = selectCuts(base, cutsPerBin, maxCuts, seed);
    vector<ScoredMove> screen;
    int totalTargets = 0;
    for (size_t ci = 0; ci < cuts.size(); ++ci) {
        SPRMove const cut = cuts[ci];
        vector<SPRMove> moves = allMovesForCut(base, cut.cutParent, cut.cutChild,
                                               maxTargets, seed + ci + 1);
        totalTargets += (int)moves.size();
        auto batch = generalizedBatch(base, cut, genes, Tscreen, threads,
                                      oldScreen, (int)nTaxa);
        int const bucket = cutBin(base.sz[cut.cutChild]);
        for (auto const& m : moves)
            screen.push_back({m, bucket, base.sz[cut.cutChild],
                              batchDelta(batch, m.targetChild)});
        std::cerr << "[spr-search] cut " << ci << " targets=" << moves.size() << '\n';
    }
    std::array<vector<ScoredMove>, 8> bins;
    for (auto const& x : screen) bins[(size_t)x.bucket].push_back(x);
    std::set<MoveKey> shortlistKeys;
    vector<ScoredMove> shortlist;
    auto addShort = [&](ScoredMove const& x) {
        if (shortlistKeys.insert(keyOf(x.move)).second) shortlist.push_back(x);
    };
    for (auto& bin : bins) {
        std::sort(bin.begin(), bin.end(), [](auto const& a, auto const& b) { return a.delta > b.delta; });
        for (int i = 0; i < std::min<int>(perBin, bin.size()); ++i) addShort(bin[(size_t)i]);
    }
    std::sort(screen.begin(), screen.end(), [](auto const& a, auto const& b) { return a.delta > b.delta; });
    for (int i = 0; i < std::min<int>(globalTop, screen.size()); ++i) addShort(screen[(size_t)i]);
    std::cerr << "[spr-search] cuts=" << cuts.size() << " targets=" << totalTargets
              << " shortlist=" << shortlist.size() << '\n';

    vector<ScoredMove> full;
    std::map<std::tuple<int,int>, BatchResult> batches;
    for (auto const& x : shortlist) {
        auto key = std::make_tuple(x.move.cutParent, x.move.cutChild);
        if (!batches.contains(key))
            batches.emplace(key, generalizedBatch(base, x.move, genes, Tfull, threads,
                                                  oldFull, (int)nTaxa));
        auto const& b = batches.at(key);
        __int128 const d = batchDelta(b, x.move.targetChild);
        if (!validateCandidate(base, x.move, genes, Tfull, threads, oldFull, d)) {
            std::cerr << "FAIL candidate oracle validation\n";
            return 1;
        }
        full.push_back({x.move, x.bucket, x.cutSize, d});
        std::cerr << "[spr-search] full bucket=" << cutBinName(x.bucket)
                  << " cut=" << x.cutSize << " delta_score=" << printInt128(d / 2)
                  << " screen_score=" << printInt128(x.delta / 2) << '\n';
    }
    auto best = std::max_element(full.begin(), full.end(),
        [](auto const& a, auto const& b) { return a.delta < b.delta; });
    __int128 bestDelta = full.empty() ? 0 : best->delta;
    if (best == full.end() || bestDelta <= 0) {
        std::cerr << "[spr-search] no strictly improving SPR\n";
        if (!outputFile.empty()) { std::ofstream f(outputFile); f << toNewick(base, base.root) << ";\n"; }
        return 0;
    }
    applySPR(base, best->move);
    if (!checkInvariants(base, (int)nTaxa)) return 1;
    std::cerr << "[spr-search] accepted delta_score=" << printInt128(bestDelta / 2)
              << " new_score=" << printInt128(sumCache(refreshCache(base, genes, Tfull, threads)) / 2) << '\n';
    if (!outputFile.empty()) { std::ofstream f(outputFile); f << toNewick(base, base.root) << ";\n"; }
    return 0;
}
