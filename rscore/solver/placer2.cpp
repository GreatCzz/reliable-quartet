// placer2: 用"扇区聚合公式"一次算出所有插入位置的 Δ（O(n·G)），逐类群贪心放置。
//
// 用法:
//   ./placer2 -i genes.trees -o out.tree [-t THREADS] [--genes K] [--order freq|rand] [--seed S]

#include "placedelta.hpp"

static unordered_map<string,int> __name2id_local;

int main(int argc, char** argv) {
    string genesFile, outFile;
    int nGenes = 20, threads = 1, backbone = 0, stopAfter = 0;
    string resumeInit;
    long seed = -1;
    bool randomOrder = false;
    for (int i = 1; i < argc; ++i) {
        string a = argv[i];
        auto nx = [&]() -> string { if (i + 1 < argc) return argv[++i]; std::exit(2); };
        if (a == "-i") genesFile = nx();
        else if (a == "-o") outFile = nx();
        else if (a == "-t") threads = atoi(nx().c_str());
        else if (a == "--genes") nGenes = atoi(nx().c_str());
        else if (a == "--seed") seed = atol(nx().c_str());
        else if (a == "--backbone") backbone = atoi(nx().c_str());
        else if (a == "--stop-after") stopAfter = atoi(nx().c_str());
        else if (a == "--resume-init") resumeInit = nx();
        else if (a == "--order") { string v = nx(); randomOrder = (v == "rand"); }
        else if (a == "-h") { cout << "usage: placer2 -i genes.trees -o out.tree [--genes K] [--order freq|rand] [--seed S]\n"; return 0; }
        else { cerr << "unknown arg " << a << "\n"; return 2; }
    }
    if (genesFile.empty()) { cerr << "need -i\n"; return 2; }

    vector<string> geneLines;
    if (!readAllLines(genesFile, geneLines)) { cerr << "cannot open " << genesFile << "\n"; return 1; }
    if (geneLines.empty()) { cerr << "no genes\n"; return 1; }
    int G = std::min<int>(nGenes, (int)geneLines.size());

    // 类群表与频次
    unordered_map<string,int> name2id;
    vector<string> names;
    vector<long long> freq;
    for (auto const& line : geneLines) {
        GTree g;
        int i = 0; while (i < (int)line.size() && line[i] != '(') ++i;
        if (i >= (int)line.size()) continue;
        g.root = parseNode(line, i, g);
        for (int u = 0; u < (int)g.ch.size(); ++u) if (g.ch[u].empty()) {
            auto it = name2id.find(g.lab[u]);
            int id;
            if (it == name2id.end()) { id = (int)names.size(); name2id.emplace(g.lab[u], id); names.push_back(g.lab[u]); freq.push_back(0); }
            else id = it->second;
            freq[(size_t)id]++;
        }
    }
    int const nTaxa = (int)names.size();
    cerr << "[placer2] taxa=" << nTaxa << " genes=" << G << "\n";

    __name2id_local = name2id;
    vector<int> order((size_t)nTaxa);
    for (int i = 0; i < nTaxa; ++i) order[(size_t)i] = i;
    if (!randomOrder) sort(order.begin(), order.end(), [&](int a, int b) { return freq[(size_t)a] > freq[(size_t)b]; });
    if (seed >= 0) { std::mt19937_64 rng((uint64_t)seed); std::shuffle(order.begin(), order.end(), rng); }

    PTree T;
    T.name = names;
    vector<char> placed((size_t)nTaxa, 0);
    auto addLeaf = [&](int x) { int id = T.add(-1, -1, -1, x); placed[(size_t)x] = 1; return id; };
    if (!resumeInit.empty()) {
        vector<string> lines; readAllLines(resumeInit, lines);
        if (lines.empty()) { cerr << "cannot read " << resumeInit << "\n"; return 1; }
        GTree g; int i = 0; while (i < (int)lines[0].size() && lines[0][i] != '(') ++i;
        if (i < (int)lines[0].size()) g.root = parseNode(lines[0], i, g);
        // 转换 GTree -> PTree（按名字查 id）

        function<int(int)> conv = [&](int u) -> int {
            if (g.ch[u].empty()) {
                auto it = __name2id_local.find(g.lab[u]);
                if (it == __name2id_local.end()) { cerr << "unknown taxon in init: " << g.lab[u] << "\n"; exit(1); }
                return addLeaf(it->second);
            }
            int a = conv(g.ch[u][0]), b = conv(g.ch[u][1]);
            int id = T.add(-1, a, b, -1); T.par[a] = id; T.par[b] = id; return id;
        };
        int rootId = conv(g.root);
        T.root = rootId;
        cerr << "[placer2] resume from " << resumeInit << ", placed " << (int)T.par.size() << " nodes\n";
    }
    PTree frozen;            // 骨架快照（两段式放置用）
    bool haveFrozen = false;
    auto t0 = chrono::steady_clock::now();

    int nPlaced = (int)count(placed.begin(), placed.end(), (char)1);
    for (int idx = 0; idx < nTaxa; ++idx) {
        int const x = order[(size_t)idx];
        if (placed[(size_t)x]) continue;
        if (stopAfter > 0 && nPlaced >= stopAfter) break;
        nPlaced++;
        if (T.par.empty()) { T.add(-1, -1, -1, x); T.root = 0; placed[(size_t)x] = 1; continue; }
        if ((int)T.par.size() == 1) {
            int p = T.root;
            int leaf = T.add(T.par[p], -1, -1, x);
            int w = T.add(T.par[p], p, leaf, -1);
            T.par[p] = w; T.par[leaf] = w; T.root = w; placed[(size_t)x] = 1;
            continue;
        }
        if (backbone > 0 && idx >= backbone && !haveFrozen) { frozen = T; haveFrozen = true; }
        PTree const& evalTree = (backbone > 0 && idx >= backbone) ? frozen : T;
        vector<long long> delta; long long rootDelta = 0;
        computeAllDeltas(evalTree, geneLines, G, names[(size_t)x], delta, rootDelta);

        long long best = rootDelta; int bestP = -1;   // -1 表示"插到根之上"
        for (int u = 0; u < (int)evalTree.par.size(); ++u) {
            if (u == evalTree.root) continue;
            if (delta[(size_t)u] > best) { best = delta[(size_t)u]; bestP = u; }
        }
        if (bestP < 0) {   // 根之上：新根 = (旧根, x)
            int leaf = T.add(-1, -1, -1, x);
            int w = T.add(-1, T.root, leaf, -1);
            T.par[leaf] = w; T.par[T.root] = w; T.root = w;
        } else {
            int p = bestP;
            int par0 = T.par[p];
            int leaf = T.add(par0, -1, -1, x);
            int w = T.add(par0, p, leaf, -1);
            if (par0 < 0) T.root = w;
            else { if (T.c0[par0] == p) T.c0[par0] = w; else T.c1[par0] = w; }
            T.par[p] = w;
            T.par[leaf] = w;
        }
        if ((idx + 1) % 200 == 0 || idx + 1 == nTaxa)
            cerr << "[placer2] " << (idx + 1) << "/" << nTaxa << " ("
                 << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - t0).count() << " s)\n";
    }

    function<string(int)> emit = [&](int u) -> string {
        if (T.taxon[u] >= 0) return names[(size_t)T.taxon[u]];
        return "(" + emit(T.c0[u]) + "," + emit(T.c1[u]) + ")";
    };
    string nw = emit(T.root) + ";";
    if (outFile.empty()) cout << nw << "\n";
    else { ofstream f(outFile); f << nw << "\n"; }
    cerr << "[placer2] done\n";
    return 0;
}
