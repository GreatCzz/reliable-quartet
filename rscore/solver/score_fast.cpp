// solution_fast.cpp
// C++ port of solution_fast.py (+ count_quartets_all.py 的预处理流程)
// 公共四叶树计数 (D/S/C 版), O(nH log n)
// 用法: ./solution_fast <species.trees> <gene.trees> [LIMIT] [THREADS]
//   物种树取文件第一非空行; 基因树逐行处理, LIMIT>0 时只跑前 LIMIT 棵(默认全部).
//   多叉 resolve 默认用"按子树叶数为权重的哈夫曼"; 环境变量 RESOLVE_MODE=balanced|huffman|sequential 可覆盖.
//   STATS=1 输出叶节点深度统计; DBGLEAF=<叶名> 打印该叶深度.
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;

// 每棵树的答案 ≤ C(n,4) 量级, int64 够用; 但 1000 棵求和需 128 位累加器
static std::string i128str(__int128 v) {
    if (v == 0) return "0";
    bool neg = v < 0;
    unsigned __int128 u = neg ? (unsigned __int128)(-v) : (unsigned __int128)v;
    char buf[48]; int p = 48;
    while (u) { buf[--p] = char('0' + (int)(u % 10)); u /= 10; }
    if (neg) buf[--p] = '-';
    return std::string(buf + p, buf + 48);
}

// ---------------- 节点池 ----------------
struct Node { vector<int> ch; string label; };
static thread_local vector<Node> pool;   // 每线程一份

static inline int newNode() { pool.push_back(Node()); return (int)pool.size() - 1; }
static inline int codeOf(const string& s) {
    if (s == "C") return 1; if (s == "D") return 2; if (s == "S") return 3; return 0;
}

// 去枝长 ":...."
static string stripBranches(const string& s) {
    string out; out.reserve(s.size());
    size_t i = 0, n = s.size();
    while (i < n) {
        char c = s[i];
        if (c == ':') { ++i; while (i < n && s[i] != ',' && s[i] != '(' && s[i] != ')' && s[i] != ';') ++i; }
        else { out.push_back(c); ++i; }
    }
    return out;
}

// 迭代 Newick 解析(支持多叉), 返回根 id
static int parse(const string& s0) {
    string s = stripBranches(s0);
    int n = (int)s.size(), i = 0;
    vector<int> st; int root = -1;
    auto isDelim = [](char c){ return c==',' || c=='(' || c==')' || c==';'; };
    while (i < n) {
        char c = s[i];
        if (c == '(') {
            int id = newNode();
            if (!st.empty()) pool[st.back()].ch.push_back(id); else root = id;
            st.push_back(id); ++i;
        } else if (c == ',') {
            ++i;
        } else if (c == ')') {
            ++i; int j = i; while (j < n && !isDelim(s[j])) ++j;
            if (!st.empty()) { pool[st.back()].label = s.substr(i, j - i); st.pop_back(); }
            i = j;
        } else if (c == ';') {
            ++i;
        } else {
            int j = i; while (j < n && !isDelim(s[j])) ++j;
            int id = newNode(); pool[id].label = s.substr(i, j - i);
            if (!st.empty()) pool[st.back()].ch.push_back(id); else root = id;
            i = j;
        }
    }
    return root;
}

// 通用后序(节点 id)
static vector<int> postorderNodes(int root) {
    vector<int> order, st{root};
    while (!st.empty()) { int t = st.back(); st.pop_back(); order.push_back(t); for (int c : pool[t].ch) st.push_back(c); }
    reverse(order.begin(), order.end());
    return order;
}

// 以 target 叶的邻接节点为根重建(返回新根 id); target 不存在则原样返回
static bool g_lastRerootFound = false;
static int rerootAtLeaf(int root, const string& target) {
    vector<vector<int>> adj(pool.size());
    vector<vector<string>> adjLab(pool.size());
    unordered_map<string,int> lab2id;
    vector<pair<int,int>> stk; stk.push_back({root, -1});
    while (!stk.empty()) {
        auto [nd, pid] = stk.back(); stk.pop_back();
        if (pid != -1) {
            adj[nd].push_back(pid); adjLab[nd].push_back(pool[nd].label);
            adj[pid].push_back(nd); adjLab[pid].push_back(pool[nd].label);
        }
        if (pool[nd].ch.empty()) lab2id[pool[nd].label] = nd;
        for (int c : pool[nd].ch) stk.push_back({c, nd});
    }
    auto it = lab2id.find(target);
    if (it == lab2id.end()) { g_lastRerootFound = false; return root; }
    g_lastRerootFound = true;
    int tgt = it->second, nroot = adj[tgt][0];
    int nrootNew = newNode(); pool[nrootNew].label = "";
    vector<tuple<int,int,int>> st2; st2.push_back({nroot, -1, nrootNew});
    while (!st2.empty()) {
        auto [cur, p, nid] = st2.back(); st2.pop_back();
        for (size_t k = 0; k < adj[cur].size(); ++k) {
            int nb = adj[cur][k]; if (nb == p) continue;
            int cid = newNode(); pool[cid].label = adjLab[cur][k];
            pool[nid].ch.push_back(cid);
            st2.push_back({nb, cur, cid});
        }
    }
    return nrootNew;
}

static void renameLeaves(int root, const function<string(const string&)>& fn) {
    auto order = postorderNodes(root);
    for (int t : order) if (pool[t].ch.empty()) pool[t].label = fn(pool[t].label);
}

static inline int pri(const string& l) {
    if (l == "D") return 3; if (l == "S") return 2; if (l == "C") return 1; return 2;
}
static void splice1(int root) {
    auto order = postorderNodes(root);
    for (int t : order) {
        vector<int> newch; newch.reserve(pool[t].ch.size());
        for (int c : pool[t].ch) {
            while (pool[c].ch.size() == 1) {
                int g = pool[c].ch[0];
                if (pri(pool[c].label) >= pri(pool[g].label)) pool[g].label = pool[c].label;
                c = g;
            }
            newch.push_back(c);
        }
        pool[t].ch.swap(newch);
    }
}

static int mkBalanced(const vector<int>& kids, int a, int b) {
    if (b - a == 1) return kids[a];
    int m = (a + b) / 2;
    int id = newNode(); pool[id].label = "C";
    int L = mkBalanced(kids, a, m), R = mkBalanced(kids, m, b);
    pool[id].ch.push_back(L); pool[id].ch.push_back(R);
    return id;
}
static void resolveBalanced(int root) {
    auto order = postorderNodes(root);
    for (int t : order) {
        if (pool[t].ch.size() > 2) {
            vector<int> kids = pool[t].ch;
            int k = (int)kids.size(), h = k / 2;
            int L = mkBalanced(kids, 0, h), R = mkBalanced(kids, h, k);
            pool[t].ch.clear(); pool[t].ch.push_back(L); pool[t].ch.push_back(R);
        }
    }
}

// 按"子树叶数"为权重的哈夫曼 resolve: 每次合并叶数最小的两个
static void resolveHuffman(int root) {
    auto order = postorderNodes(root);
    vector<int> sz(pool.size(), 1);
    for (int t : order) {
        auto& ch = pool[t].ch;
        if (ch.empty()) { sz[t] = 1; continue; }
        int s = 0; for (int c : ch) s += sz[c];
        sz[t] = s;
    }
    for (int t : order) {
        auto& ch = pool[t].ch;
        if (ch.size() <= 2) continue;
        priority_queue<pair<int,int>, vector<pair<int,int>>, greater<>> pq;
        for (int c : ch) pq.push({sz[c], c});
        while (pq.size() > 2) {
            auto [sa, a] = pq.top(); pq.pop();
            auto [sb, b] = pq.top(); pq.pop();
            int id = newNode(); pool[id].label = "C";
            pool[id].ch.push_back(a); pool[id].ch.push_back(b);
            sz.push_back(sa + sb);
            pq.push({sa + sb, id});
        }
        auto [s1, a] = pq.top(); pq.pop();
        auto [s2, b] = pq.top(); pq.pop();
        pool[t].ch.clear(); pool[t].ch.push_back(a); pool[t].ch.push_back(b);
    }
}

// 顺序链式 resolve: (((G1,G2),G3),G4)  —— 用于对拍"形状无关性"
static void resolveSequential(int root) {
    auto order = postorderNodes(root);
    for (int t : order) {
        auto& ch = pool[t].ch;
        if (ch.size() <= 2) continue;
        vector<int> kids = ch;
        int cur = kids[0];
        for (int i = 1; i + 1 < (int)kids.size(); ++i) {
            int id = newNode(); pool[id].label = "C";
            pool[id].ch.push_back(cur); pool[id].ch.push_back(kids[i]);
            cur = id;
        }
        pool[t].ch.clear(); pool[t].ch.push_back(cur); pool[t].ch.push_back(kids.back());
    }
}

static int g_resolveMode = 1;   // 0=balanced 1=huffman(叶数, 默认) 2=sequential
static bool g_stats = false;    // 统计叶节点平均深度
static atomic<long long> g_depthSum{0}, g_leafN{0}, g_depthMin{LLONG_MAX}, g_depthMax{0};

// ---------------- 扁平化 ----------------
struct Flat {
    vector<int> lc, rc, par, code;
    vector<char> isleaf;
    unordered_map<string,int> leafid;
    int rootid = -1;
};
static Flat flatten(int root) {
    Flat F;
    auto newn = [&](int code)->int {
        F.lc.push_back(-1); F.rc.push_back(-1); F.par.push_back(-1);
        F.code.push_back(code); F.isleaf.push_back(0);
        return (int)F.lc.size() - 1;
    };
    F.rootid = newn(codeOf(pool[root].label));
    vector<pair<int,int>> stk; stk.push_back({root, F.rootid});
    while (!stk.empty()) {
        auto [nd, i] = stk.back(); stk.pop_back();
        if (!pool[nd].ch.empty()) {
            vector<int> ids;
            for (int c : pool[nd].ch) { int ci = newn(codeOf(pool[c].label)); F.par[ci] = i; ids.push_back(ci); }
            F.lc[i] = ids[0]; F.rc[i] = ids[1];
            for (size_t k = 0; k < pool[nd].ch.size(); ++k) stk.push_back({pool[nd].ch[k], ids[k]});
        } else { F.isleaf[i] = 1; F.leafid[pool[nd].label] = i; }
    }
    return F;
}
static vector<int> flatPostorder(const vector<int>& lc, const vector<int>& rc, int rootid) {
    vector<int> order, st{rootid};
    while (!st.empty()) { int t = st.back(); st.pop_back(); order.push_back(t); if (lc[t] >= 0) { st.push_back(lc[t]); st.push_back(rc[t]); } }
    reverse(order.begin(), order.end());
    return order;
}

// ---------------- 全局 DP 状态 (每线程一份) ----------------
static thread_local vector<array<ll,13>> gV1, gV2;
static thread_local vector<char> gcu, gnu;
static thread_local vector<int> glc, grc, gpar;
static thread_local int ggeRoot;
static thread_local vector<vector<int>> gspCh, gcache;

static inline void step1(int t) {
    int l = glc[t], r = grc[t];
    ll nu = gnu[l], nv = gnu[r];
    auto& U = gV1[l]; auto& V = gV1[r]; auto& A = gV1[t];
    ll Ru=U[0],Bu=U[1],Gu=U[2],RRu=U[3],BBu=U[4],GGu=U[5],RBu=U[6],RGu=U[7],BGu=U[8],coRu=U[9],coBu=U[10],coGu=U[11],Qu=U[12];
    ll Rv=V[0],Bv=V[1],Gv=V[2],RRv=V[3],BBv=V[4],GGv=V[5],RBv=V[6],RGv=V[7],BGv=V[8],coRv=V[9],coBv=V[10],coGv=V[11],Qv=V[12];
    A[0]=nu*Ru+nv*Rv; A[1]=nu*Bu+nv*Bv; A[2]=nu*Gu+nv*Gv;
    A[3]=RRu+RRv+nu*nv*Ru*Rv; A[4]=BBu+BBv+nu*nv*Bu*Bv; A[5]=GGu+GGv+nu*nv*Gu*Gv;
    A[6]=RBu+RBv+nu*nv*(Ru*Bv+Bu*Rv);
    A[7]=RGu+RGv+nu*nv*(Ru*Gv+Gu*Rv);
    A[8]=BGu+BGv+nu*nv*(Bu*Gv+Gu*Bv);
    A[9]=nu*coRu+nv*coRv + (BBu*nv*Gv + BBv*nu*Gu + GGu*nv*Bv + GGv*nu*Bu + BGu*nv*Rv + BGv*nu*Ru);
    A[10]=nu*coBu+nv*coBv + (RRu*nv*Gv + RRv*nu*Gu + GGu*nv*Rv + GGv*nu*Ru + RGu*nv*Bv + RGv*nu*Bu);
    A[11]=nu*coGu+nv*coGv + (BBu*nv*Rv + BBv*nu*Ru + RRu*nv*Bv + RRv*nu*Bu + RBu*nv*Gv + RBv*nu*Gu);
    A[12]=Qu+Qv + nu*nv*(coRu*Rv + coRv*Ru + coBu*Bv + coBv*Bu + coGu*Gv + coGv*Gu)
          + (RRu*BGv + RRv*BGu + BBu*RGv + BBv*RGu + GGu*RBv + GGv*RBu);
}

static inline void step2(int t) {
    int l = glc[t], r = grc[t];
    ll cu = gcu[l], nu = gnu[l], cv = gcu[r], nv = gnu[r];
    auto& U = gV2[l]; auto& V = gV2[r]; auto& B = gV2[t];
    ll Ru=U[0],Bu=U[1],Gu=U[2],RRu=U[3],BBu=U[4],GGu=U[5],RBu=U[6],RGu=U[7],BGu=U[8],coRu=U[9],coBu=U[10],coGu=U[11],Qu=U[12];
    ll Rv=V[0],Bv=V[1],Gv=V[2],RRv=V[3],BBv=V[4],GGv=V[5],RBv=V[6],RGv=V[7],BGv=V[8],coRv=V[9],coBv=V[10],coGv=V[11],Qv=V[12];
    B[0]=nu*Ru+nv*Rv; B[1]=nu*Bu+nv*Bv; B[2]=nu*Gu+nv*Gv;
    B[3]=cu*RRu+cv*RRv+nu*nv*Ru*Rv; B[4]=cu*BBu+cv*BBv+nu*nv*Bu*Bv; B[5]=cu*GGu+cv*GGv+nu*nv*Gu*Gv;
    B[6]=cu*RBu+cv*RBv+nu*nv*(Ru*Bv+Bu*Rv);
    B[7]=cu*RGu+cv*RGv+nu*nv*(Ru*Gv+Gu*Rv);
    B[8]=cu*BGu+cv*BGv+nu*nv*(Bu*Gv+Gu*Bv);
    B[9]=nu*coRu+nv*coRv + cu*nv*(BBu*Gv + GGu*Bv + BGu*Rv) + cv*nu*(BBv*Gu + GGv*Bu + BGv*Ru);
    B[10]=nu*coBu+nv*coBv + cu*nv*(RRu*Gv + GGu*Rv + RGu*Bv) + cv*nu*(RRv*Gu + GGv*Ru + RGv*Bu);
    B[11]=nu*coGu+nv*coGv + cu*nv*(RRu*Bv + BBu*Rv + RBu*Gv) + cv*nu*(RRv*Bu + BBv*Ru + RBv*Gu);
    B[12]=Qu+Qv + nu*nv*(coRu*Rv + coRv*Ru + coBu*Bv + coBv*Bu + coGu*Gv + coGv*Gu)
          + cu*cv*(RRu*BGv + RRv*BGu + BBu*RGv + BBv*RGu + GGu*RBv + GGv*RBu);
}

static void colorLeaf(int v, int c) {  // c: 0=R,1=B,2=G
    auto& A = gV1[v]; auto& B = gV2[v];
    A[0] = (c==0); A[1] = (c==1); A[2] = (c==2);
    B[0] = (c==0); B[1] = (c==1); B[2] = (c==2);
    int t = gpar[v];
    while (t >= 0) { step1(t); step2(t); t = gpar[t]; }
}

static pair<ll,ll> cq(int R) {
    int LC = gspCh[R][0], RC = gspCh[R][1];
    const vector<int>& leafR = gcache[RC];
    ll s1 = 0, s2 = 0;
    if (!gspCh[RC].empty()) {
        auto pr = cq(RC); s1 += pr.first; s2 += pr.second;
        for (int e : leafR) colorLeaf(e, 2);
    }
    if (!gspCh[LC].empty()) {
        auto pr = cq(LC); s1 += pr.first; s2 += pr.second;
        for (int e : gcache[gspCh[LC][1]]) colorLeaf(e, 1);
    } else if (!gcache[LC].empty()) {      // LC 是叶子；该基因树可能缺这个类群
        colorLeaf(gcache[LC][0], 1);
    }
    for (int e : leafR) colorLeaf(e, 0);
    return {s1 + gV1[ggeRoot][12], s2 + gV2[ggeRoot][12]};
}

// 物种树: 后序计算叶数并 left-heavier 重排 children, 同时建 cache(sp节点 -> 基因叶id)
static void spBuildCacheAndReorder(int spRoot, const Flat& F) {
    int P = (int)pool.size();
    auto order = postorderNodes(spRoot);
    vector<int> sz(P, 1);
    gcache.assign(P, {});
    for (int t : order) {
        auto& ch = pool[t].ch;
        if (ch.empty()) {
            sz[t] = 1;
            auto itLeaf = F.leafid.find(pool[t].label);
            if (itLeaf == F.leafid.end()) gcache[t] = {};      // 该基因树缺这个类群
            else gcache[t] = { itLeaf->second };
            continue;
        }
        vector<int> list;
        for (int c : ch) for (int x : gcache[c]) list.push_back(x);
        int a = ch[0], b = ch[1];
        if (sz[a] < sz[b]) { swap(ch[0], ch[1]); swap(a, b); }
        sz[t] = sz[ch[0]] + sz[ch[1]];
        gcache[t] = std::move(list);
    }
    gspCh.assign(P, {});
    for (int t : order) if (!pool[t].ch.empty()) gspCh[t] = pool[t].ch;
}

// 物种树只读共享副本
static vector<Node> spPool;
static int spRootBase = -1;
static int spN = 0;

static thread_local int g_lastLeafCount = 0;
static std::unordered_set<std::string> g_geneTaxa;
static std::mutex g_taxaMutex;

// 处理一棵基因树(线程内), 返回其公共四叶树数
static ll solveGene(const string& gline) {
    if (pool.size() < (size_t)spN) pool = spPool; else pool.resize(spN);  // 复用物种树节点
    int geRaw = parse(gline);
    geRaw = rerootAtLeaf(geRaw, "1_0_0");
    renameLeaves(geRaw, [](const string& nm){ size_t p = nm.find('_'); return nm.substr(0, p); });
    splice1(geRaw);
    while ((int)pool[geRaw].ch.size() == 1) {      // 根单叉：直接提升唯一孩子
        int const g = pool[geRaw].ch[0];
        pool[g].label = pool[geRaw].label;
        geRaw = g;
    }
    if (g_resolveMode == 1) resolveHuffman(geRaw);
    else if (g_resolveMode == 2) resolveSequential(geRaw);
    else resolveBalanced(geRaw);
    Flat F = flatten(geRaw);

    vector<int> order = flatPostorder(F.lc, F.rc, F.rootid);
    {   // gene left-heavy
        int N = (int)F.lc.size(); vector<int> sz(N, 0);
        for (int t : order) {
            if (F.lc[t] < 0) sz[t] = 1;
            else {
                int l = F.lc[t], r = F.rc[t];
                if (sz[l] < sz[r]) { swap(F.lc[t], F.rc[t]); l = F.lc[t]; r = F.rc[t]; }
                sz[t] = sz[l] + sz[r];
            }
        }
    }

    if (g_stats) {   // 统计叶节点平均深度(resolve 后, 根深度=0)
        int NN = (int)F.lc.size();
        vector<int> dep(NN, 0);
        long long dsum = 0, ln = 0, dmin = LLONG_MAX, dmax = 0;
        for (int i = NN - 1; i >= 0; --i) {       // order 为后序, 反向 = 自顶向下
            int t = order[i];
            if (F.lc[t] < 0) { dsum += dep[t]; ++ln; if (dep[t] < dmin) dmin = dep[t]; if (dep[t] > dmax) dmax = dep[t]; }
            else { dep[F.lc[t]] = dep[t] + 1; dep[F.rc[t]] = dep[t] + 1; }
        }
        g_depthSum.fetch_add(dsum, memory_order_relaxed);
        g_leafN.fetch_add(ln, memory_order_relaxed);
        long long cur = g_depthMin.load();
        while (dmin < cur && !g_depthMin.compare_exchange_weak(cur, dmin)) {}
        cur = g_depthMax.load();
        while (dmax > cur && !g_depthMax.compare_exchange_weak(cur, dmax)) {}
        if (const char* dbg = getenv("DBGLEAF")) {
            auto it = F.leafid.find(dbg);
            fprintf(stderr, "  [dbg] leaf %s depth=%d (rerootFound=%d, rootid=%d, N=%d, orderN=%zu)\n",
                    dbg, it != F.leafid.end() ? dep[it->second] : -1, (int)g_lastRerootFound, F.rootid, NN, order.size());
        }
    }

    int N = (int)F.lc.size();
    glc = F.lc; grc = F.rc; gpar = F.par; ggeRoot = F.rootid;
    gcu.assign(N, 0); gnu.assign(N, 0);
    for (int t = 0; t < N; ++t) { gcu[t] = (F.code[t] == 1); gnu[t] = (F.code[t] != 2); }
    gV1.assign(N, array<ll,13>{}); gV2.assign(N, array<ll,13>{});
    for (int t = 0; t < N; ++t) if (F.isleaf[t]) { gV1[t][2] = 1; gV2[t][2] = 1; }
    for (int t : order) if (glc[t] >= 0) { step1(t); step2(t); }

    g_lastLeafCount = (int)F.leafid.size();
    { std::lock_guard<std::mutex> lk(g_taxaMutex);
      for (auto const& kv : F.leafid) g_geneTaxa.insert(kv.first); }

    spBuildCacheAndReorder(spRootBase, F);
    auto pr = cq(spRootBase);
    return pr.first / 2 - pr.second / 2;
}

static bool readAllLines(const string& path, vector<string>& out) {
    ifstream f(path);
    if (!f.is_open()) return false;
    string l;
    while (getline(f, l)) {
        while (!l.empty() && (l.back() == '\r' || l.back() == ' ' || l.back() == '\t')) l.pop_back();
        bool nz = false;
        for (char c : l) if (!isspace((unsigned char)c)) { nz = true; break; }
        if (nz) out.push_back(l);
    }
    return true;
}

int main(int argc, char** argv) {
    string treeFile, genesFile;
    int threads = 0, limit = 0;
    bool quiet = false;
    for (int i = 1; i < argc; ++i) {
        string a = argv[i];
        auto nx = [&]() -> string { if (i + 1 < argc) return string(argv[++i]); std::exit(2); };
        if (a == "-t" || a == "--tree") treeFile = nx();
        else if (a == "-g" || a == "--genes") genesFile = nx();
        else if (a == "-j" || a == "--threads") threads = atoi(nx().c_str());
        else if (a == "-q" || a == "--quiet") quiet = true;
        else if (a == "--limit") limit = atoi(nx().c_str());
        else if (a == "-h" || a == "--help") {
            cout << "usage: score_fast -t <species.tree> -g <genes.trees> [-j THREADS] [-q] [--limit N]\n";
            return 0;
        }
        else { cerr << "Error: unknown argument " << a << "\n"; return 2; }
    }
    if (treeFile.empty() || genesFile.empty()) {
        cerr << "Error: -t <species tree> and -g <gene trees> are required\n";
        return 2;
    }
    if (const char* rm = getenv("RESOLVE_MODE")) {
        if (!strcmp(rm, "huffman")) g_resolveMode = 1;
        else if (!strcmp(rm, "sequential")) g_resolveMode = 2;
        else g_resolveMode = 0;
    }

    vector<string> spLines, geneLines;
    if (!readAllLines(treeFile, spLines)) { cerr << "Error: cannot open species tree file: " << treeFile << "\n"; return 1; }
    if (spLines.empty()) { cerr << "Error: no species tree found in " << treeFile << "\n"; return 1; }
    if (spLines.size() != 1) {
        cerr << "Error: expected exactly one species tree in " << treeFile << ", found " << spLines.size() << "\n";
        return 1;
    }
    if (!readAllLines(genesFile, geneLines)) { cerr << "Error: cannot open gene tree file: " << genesFile << "\n"; return 1; }
    if (geneLines.empty()) { cerr << "Error: no gene tree found in " << genesFile << "\n"; return 1; }
    if (limit > 0 && (int)geneLines.size() > limit) geneLines.resize(limit);

    int r = 0;
    try { r = parse(spLines[0]); }
    catch (std::exception const& e) { cerr << "Error: bad Newick in species tree: " << e.what() << "\n"; return 1; }
    if (pool[r].ch.size() == 3) {           // 根三叉 = 无根二叉树的标准写法
        int synth = newNode();
        pool[synth].label = "";
        pool[synth].ch.push_back(pool[r].ch[1]);
        pool[synth].ch.push_back(pool[r].ch[2]);
        pool[r].ch.resize(1);
        pool[r].ch.push_back(synth);
    } else if (pool[r].ch.size() != 2 && !pool[r].ch.empty()) {
        cerr << "Error: species tree must be binary (root with " << pool[r].ch.size() << " children)\n";
        return 1;
    }
    {   // 检查所有内部节点二叉
        vector<int> st{r};
        while (!st.empty()) {
            int t = st.back(); st.pop_back();
            if (pool[t].ch.empty()) continue;
            if (pool[t].ch.size() != 2) {
                cerr << "Error: species tree must be binary (node with " << pool[t].ch.size() << " children)\n";
                return 1;
            }
            for (int c : pool[t].ch) st.push_back(c);
        }
    }
    spRootBase = r;
    spN = (int)pool.size();
    spPool = pool;

    // 物种树叶名集合与叶数
    std::unordered_set<std::string> spTaxa;
    {
        vector<int> st{r};
        while (!st.empty()) {
            int t = st.back(); st.pop_back();
            if (pool[t].ch.empty()) spTaxa.insert(pool[t].label);
            else for (int c : pool[t].ch) st.push_back(c);
        }
    }
    int const spLeafCount = (int)spTaxa.size();

    int n = (int)geneLines.size();
    int nthreads = threads;
    if (nthreads <= 0) nthreads = (int)thread::hardware_concurrency();
    if (nthreads < 1) nthreads = 1;
    if (nthreads > n) nthreads = n;

    vector<ll> res(n, 0);
    vector<char> skippedFlag(n, 0);
    auto t0 = chrono::steady_clock::now();
    atomic<int> idx{0};
    auto worker = [&]() {
        int i;
        while ((i = idx.fetch_add(1, memory_order_relaxed)) < n) {
            res[i] = solveGene(geneLines[i]);
            if (g_lastLeafCount < 4) { res[i] = 0; skippedFlag[i] = 1; }
        }
    };
    vector<thread> ths;
    ths.reserve(nthreads);
    for (int t = 0; t < nthreads; ++t) ths.emplace_back(worker);
    for (auto& th : ths) th.join();
    double const secs = chrono::duration<double>(chrono::steady_clock::now() - t0).count();

    __int128 total = 0;
    int nSkipped = 0;
    for (int i = 0; i < n; ++i) { total += res[i]; if (skippedFlag[i]) ++nSkipped; }

    int nUnused = 0, nIgnored = 0;
    for (auto const& name : spTaxa) if (!g_geneTaxa.count(name)) ++nUnused;
    for (auto const& name : g_geneTaxa) if (!spTaxa.count(name)) ++nIgnored;
    if (!quiet) {
        if (nUnused) cerr << "Warning: " << nUnused << " species-tree taxa never occur in any gene tree (they contribute nothing)\n";
        if (nIgnored) cerr << "Warning: " << nIgnored << " gene-tree taxa are absent from the species tree (ignored)\n";
        if (nSkipped) cerr << "Warning: skipped " << nSkipped << " gene tree(s) with <4 taxa\n";
    }

    cout << "{\"score\":" << i128str(total)
         << ",\"n_taxa\":" << spLeafCount
         << ",\"n_genes\":" << n
         << ",\"n_scored_genes\":" << (n - nSkipped)
         << ",\"n_skipped_genes\":" << nSkipped
         << ",\"n_unused_taxa\":" << nUnused
         << ",\"n_ignored_taxa\":" << nIgnored
         << ",\"runtime_sec\":" << fixed << setprecision(3) << secs
         << ",\"per_gene\":[";
    for (int i = 0; i < n; ++i) { if (i) cout << ','; cout << res[i]; }
    cout << "]}\n";
    return 0;
}
