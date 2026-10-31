#include <array>
#include <cstdint>
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

struct BruteGene {
    vector<int> parent, leafNode;
    vector<char> parentSolid;
    vector<uint64_t> subtreeMask;
};

BruteGene prepareBruteGene(GeneTree const& g, int nTaxa) {
    BruteGene b;
    b.parent.assign(g.n, -1);
    b.parentSolid.assign(g.n, 1);
    b.leafNode.assign(nTaxa, -1);
    b.subtreeMask.assign(g.n, 0);
    for (int u = 0; u < g.n; ++u) {
        if (g.taxon[u] >= 0) b.leafNode[g.taxon[u]] = u;
        for (int k = g.start[u]; k < g.start[u + 1]; ++k) {
            int const child = g.kids[k];
            b.parent[child] = u;
            b.parentSolid[child] = g.solid[child];
        }
    }
    for (int u = g.n - 1; u >= 0; --u) {
        if (g.taxon[u] >= 0) b.subtreeMask[u] = 1ULL << g.taxon[u];
        else for (int k = g.start[u]; k < g.start[u + 1]; ++k)
            b.subtreeMask[u] |= b.subtreeMask[g.kids[k]];
    }
    return b;
}

bool bruteSolidPath(BruteGene const& g, int a, int b) {
    int u = g.leafNode[a], v = g.leafNode[b];
    if (u < 0 || v < 0) return false;
    vector<char> ancestors(g.parent.size(), 0);
    for (int x = u; x >= 0; x = g.parent[x]) ancestors[x] = 1;
    int lca = v;
    while (lca >= 0 && !ancestors[lca]) lca = g.parent[lca];
    for (int x = u; x != lca; x = g.parent[x])
        if (!g.parentSolid[x]) return false;
    for (int x = v; x != lca; x = g.parent[x])
        if (!g.parentSolid[x]) return false;
    return true;
}

bool bruteDisplays(vector<uint64_t> const& masks, vector<int> const& parent,
    int a, int b, int c, int d) {
    for (size_t u = 0; u < masks.size(); ++u) {
        if (parent[u] < 0) continue;
        uint64_t const m = masks[u];
        bool const splitAB = ((m >> a) & 1) && ((m >> b) & 1)
            && !((m >> c) & 1) && !((m >> d) & 1);
        bool const splitCD = ((m >> c) & 1) && ((m >> d) & 1)
            && !((m >> a) & 1) && !((m >> b) & 1);
        if (splitAB || splitCD) return true;
    }
    return false;
}

long long bruteOfficialScore(DynTree const& tree, vector<BruteGene> const& genes,
    int nTaxa) {
    vector<int> parent(tree.n, -1);
    vector<uint64_t> masks(tree.n, 0);
    function<uint64_t(int)> buildMask = [&](int u) -> uint64_t {
        if (tree.taxon[u] >= 0) return masks[u] = 1ULL << tree.taxon[u];
        masks[u] = buildMask(tree.c0[u]) | buildMask(tree.c1[u]);
        return masks[u];
    };
    buildMask(tree.root);
    for (int u = 0; u < tree.n; ++u) {
        if (tree.c0[u] >= 0) {
            parent[tree.c0[u]] = u;
            parent[tree.c1[u]] = u;
        }
    }

    long long score = 0;
    for (int a = 0; a < nTaxa; ++a)
    for (int b = a + 1; b < nTaxa; ++b)
    for (int c = b + 1; c < nTaxa; ++c)
    for (int d = c + 1; d < nTaxa; ++d) {
        int const splits[3][4] = {{a, b, c, d}, {a, c, b, d}, {a, d, b, c}};
        if (!((masks[tree.root] >> a) & 1) || !((masks[tree.root] >> b) & 1)
            || !((masks[tree.root] >> c) & 1) || !((masks[tree.root] >> d) & 1)) continue;
        for (auto const& gene : genes) {
            for (auto const& q : splits) {
                if (!bruteSolidPath(gene, q[0], q[1])
                    || !bruteSolidPath(gene, q[2], q[3])) continue;
                if (bruteDisplays(gene.subtreeMask, gene.parent, q[0], q[1], q[2], q[3])
                    && bruteDisplays(masks, parent, q[0], q[1], q[2], q[3])) ++score;
            }
        }
    }
    return score;
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
        "((t0,t9)S,(t2,t7)S,(t4,t6)D,t1);"
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
    vector<BruteGene> bruteGenes;
    for (auto const& gene : genes) bruteGenes.push_back(prepareBruteGene(gene, (int)nTaxa));
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
    if (tree.wsum[tree.root] / 2 != bruteOfficialScore(tree, bruteGenes, (int)nTaxa)) {
        cerr << "initial independent brute score mismatch\n";
        return 1;
    }
    checkBits(tree);

    int checkedCandidates = 0;
    int newOnlyMatches = 0;
    for (int step = 2; step < (int)nTaxa; ++step) {
        int const x = ids[(size_t)step];
        vector<u64> bX(W, 0), activePlusX = active;
        bX[(size_t)x >> 6] |= 1ULL << (x & 63);
        for (size_t w = 0; w < W; ++w) activePlusX[w] |= bX[w];
        BatchStats batch;
        vector<int> allCandidates;
        for (int p = 0; p < tree.n; ++p) allCandidates.push_back(p);
        vector<Scratch> batchScratch(1);
        batchScratch[0].resize(maxN, W);
        computeBatchStats(tree, bX.data(), activePlusX, allCandidates, genes,
            (int)genes.size(), batchScratch, 1, batch);
        long long const before = freshSum(tree, genes, active, scratch);
        long long const bruteBefore = bruteOfficialScore(tree, bruteGenes, (int)nTaxa);






        for (int p = 0; p < tree.n; ++p) {
            long long const oracle = candidateDelta(tree, p, bX.data(),
                activePlusX, genes, (int)genes.size(), scratch);
            long long const predicted = batch.base + batch.prefix[p] + batch.newNode[p];
            DynTree trial = tree;
            vector<u64> trialActive = active;
            insertTaxonAbove(trial, p, x);
            trialActive[(size_t)x >> 6] |= 1ULL << (x & 63);
            long long const actual = freshSum(trial, genes, trialActive, scratch) - before;
            long long const bruteActual = 2 * (bruteOfficialScore(trial, bruteGenes, (int)nTaxa) - bruteBefore);
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
            if (oracle != predicted || predicted != actual || predicted != bruteActual
                || predicted < 0 || (predicted & 1)) {
                cerr << "candidate delta mismatch at active size " << step
                     << ", p=" << p << ": batch=" << predicted
                     << " oracle=" << oracle << " formula=" << actual
                     << " brute=" << bruteActual << '\n';
                return 1;
            }
            checkBits(trial);
            ++checkedCandidates;
        }

        int const chosen = (step * 3) % tree.n;
        DynTree batchTree = tree;
        vector<u64> batchActive = active;
        commitBatchInsertion(batchTree, chosen, x, batchActive, batch);
        commitInsertion(tree, chosen, x, active, genes, (int)genes.size(), scratch);
        if (batchTree.n != tree.n || batchTree.root != tree.root
            || batchTree.par != tree.par || batchTree.c0 != tree.c0
            || batchTree.c1 != tree.c1 || batchTree.taxon != tree.taxon
            || batchTree.sz != tree.sz || batchTree.bits != tree.bits
            || batchTree.wsum != tree.wsum || batchActive != active) {
            cerr << "incremental cache differs from full-refresh cache at active size "
                 << step + 1 << '\n';
            return 1;
        }
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
         << checkedCandidates << " candidate deltas (formula/oracle/brute), "
         << "bitsets and caches matched; new-node-only matched "
         << newOnlyMatches << "/" << checkedCandidates << " candidates\n";
    return 0;
}
