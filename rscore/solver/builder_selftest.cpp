#include <functional>
#include <iostream>
#include <string>
#include <vector>

#define BUILDER_SELFTEST
#include "builder.cpp"

namespace {

using namespace std;

long long freshSum(DynTree const& tree, vector<GeneTree> const& genes,
    vector<u64> const& active, Scratch& scratch) {
    long long sum = 0;
    int const activeCount = tree.sz[tree.root];
    for (int v = 0; v < tree.n; ++v) {
        if (tree.c0[v] >= 0)
            sum += nodeWsum(tree, v, genes, (int)genes.size(), active,
                activeCount, scratch);
    }
    return sum;
}

void checkBits(DynTree const& tree) {
    vector<u64> expected((size_t)tree.n * tree.W, 0);
    function<void(int)> visit = [&](int u) {
        if (tree.taxon[u] >= 0) {
            size_t const word = (size_t)tree.taxon[u] >> 6;
            expected[(size_t)u * tree.W + word] |= 1ULL << (tree.taxon[u] & 63);
            return;
        }
        visit(tree.c0[u]);
        visit(tree.c1[u]);
        for (size_t w = 0; w < tree.W; ++w)
            expected[(size_t)u * tree.W + w] =
                expected[(size_t)tree.c0[u] * tree.W + w] |
                expected[(size_t)tree.c1[u] * tree.W + w];
    };
    visit(tree.root);
    for (int u = 0; u < tree.n; ++u) {
        for (size_t w = 0; w < tree.W; ++w) {
            u64 const got = tree.bits[(size_t)u * tree.W + w];
            u64 const want = expected[(size_t)u * tree.W + w];
            if (got != want) {
                cerr << "bitset mismatch at node " << u << " word " << w
                     << ": got=" << got << " want=" << want << '\n';
                exit(1);
            }
        }
    }
}

} // namespace

int main() {
    vector<string> lines = {
        "((t0,t1)S,(t2,t3)S,(t4,t5)D,(t6,t7)S,(t8,t9)S);",
        "((t0,t2)D,(t1,t4)S,(t3,t6)D,(t5,t8)S,t9);",
        "((t0,t9)S,(t2,t7)S,(t4,t6)D,t1,t3,t5);"
    };
    vector<GeneTree> genes;
    int maxN = 1;
    for (string const& line : lines) {
        newick::Newick nw(line);
        genes.push_back(buildGene(nw.root()));
        maxN = max(maxN, genes.back().n);
    }
    size_t const nTaxa = common::taxonName2ID.nTaxa();
    if (nTaxa != 10) {
        cerr << "unexpected taxon count: " << nTaxa << '\n';
        return 1;
    }
    size_t const W = (nTaxa + 63) / 64;
    Scratch scratch;
    scratch.resize(maxN, W);

    DynTree tree;
    tree.W = W;
    vector<u64> active(W, 0);
    vector<int> ids;
    for (int i = 0; i < 10; ++i)
        ids.push_back((int)common::taxonName2ID["t" + to_string(i)]);
    insertFirstTaxon(tree, ids[0]);
    active[(size_t)ids[0] >> 6] |= 1ULL << (ids[0] & 63);
    insertSecondTaxon(tree, ids[1]);
    active[(size_t)ids[1] >> 6] |= 1ULL << (ids[1] & 63);
    tree.wsum[tree.root] = freshSum(tree, genes, active, scratch);
    checkBits(tree);

    int checkedCandidates = 0;
    int newOnlyMatches = 0;
    for (int step = 2; step < (int)nTaxa; ++step) {
        int const x = ids[(size_t)step];
        vector<u64> bX(W, 0), activePlusX = active;
        bX[(size_t)x >> 6] |= 1ULL << (x & 63);
        for (size_t w = 0; w < W; ++w) activePlusX[w] |= bX[w];
        long long const before = freshSum(tree, genes, active, scratch);

        for (int p = 0; p < tree.n; ++p) {
            long long const predicted = candidateDelta(tree, p, bX.data(),
                activePlusX, genes, (int)genes.size(), scratch);
            DynTree trial = tree;
            vector<u64> trialActive = active;
            insertTaxonAbove(trial, p, x);
            trialActive[(size_t)x >> 6] |= 1ULL << (x & 63);
            long long const actual = freshSum(trial, genes, trialActive, scratch) - before;
            vector<u64> newC(W, 0);
            setComplement(newC.data(), activePlusX.data(),
                tree.bits.data() + (size_t)p * W, bX.data(), W);
            long long newOnly = 0;
            for (auto const& gene : genes)
                newOnly += Wnode(gene, tree.bits.data() + (size_t)p * W,
                    bX.data(), newC.data(), tree.sz[p], 1,
                    tree.sz[tree.root] + 1 - tree.sz[p] - 1, W, scratch);
            if (newOnly == actual) ++newOnlyMatches;
            if (newOnly != actual && checkedCandidates < 5)
                cerr << "new-node-only differs: p=" << p << " new=" << newOnly
                     << " actual=" << actual << '\n';
            if (predicted != actual) {
                cerr << "candidate delta mismatch at active size " << step
                     << ", p=" << p << ": predicted=" << predicted
                     << " actual=" << actual << '\n';
                return 1;
            }
            checkBits(trial);
            ++checkedCandidates;
        }

        int const chosen = (step * 3) % tree.n;
        commitInsertion(tree, chosen, x, active, genes, (int)genes.size(), scratch);
        long long const after = freshSum(tree, genes, active, scratch);
        long long cachedSum = 0;
        for (int v = 0; v < tree.n; ++v)
            if (tree.c0[v] >= 0) cachedSum += tree.wsum[v];
        if (cachedSum != after) {
            cerr << "cache mismatch after active size " << step + 1
                 << ": cached=" << cachedSum << " fresh=" << after << '\n';
            return 1;
        }
        checkBits(tree);
    }
    cout << "builder selftest passed: " << nTaxa << " taxa, "
         << checkedCandidates << " candidate deltas, bitsets and caches matched; "
         << "new-node-only delta matched " << newOnlyMatches << "/"
         << checkedCandidates << " candidates\n";
    return 0;
}
