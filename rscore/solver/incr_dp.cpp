// incr_dp: 13 状态 DP 原型（二叉 + C 边 + V1/V2）。
// 目的：验证 (a) 用"每节点三划分"求和能否复现总分；(b) 增量更新是否正确且更快。
// 用法: ./incr_dp -t species.tree -g genes.trees [-j N] [-q]
#include <bits/stdc++.h>
#include <omp.h>
using namespace std;
using ll = long long;
using u64 = unsigned long long;
using St = array<ll,13>;

struct Raw { vector<int> ch; string lab, name; };

static vector<Raw> parseRaw(const string& s0, int& root) {
    string s; for (char c : s0) if (c!='\r' && c!='\n') s.push_back(c);
    vector<Raw> P; vector<int> st; root=-1; int n=s.size(), i=0;
    auto isD=[](char c){ return c==','||c=='('||c==')'||c==';'; };
    while (i<n) {
        char c=s[i];
        if (c=='(') { int id=P.size(); P.push_back({}); if(!st.empty()) P[st.back()].ch.push_back(id); else root=id; st.push_back(id); ++i; }
        else if (c==',') ++i;
        else if (c==')') { ++i; int j=i; while(j<n && !isD(s[j])) ++j; if(!st.empty()){ P[st.back()].lab=s.substr(i,j-i); st.pop_back(); } i=j; }
        else if (c==';') ++i;
        else { int j=i; while(j<n && !isD(s[j])) ++j; int id=P.size(); P.push_back({}); P[id].name=s.substr(i,j-i);
            if(!st.empty()) P[st.back()].ch.push_back(id); else root=id; i=j; }
    }
    return P;
}

struct Bin { vector<int> lc, rc, par, sz, leaves, order; vector<unsigned char> code; vector<string> nm; int root=0, n=0; };

static Bin toBinary(const vector<Raw>& P, int root) {
    Bin b;
    function<int(int)> conv = [&](int u)->int {
        int id = b.n++;
        b.lc.push_back(-1); b.rc.push_back(-1); b.par.push_back(-1); b.sz.push_back(0);
        b.code.push_back(P[u].ch.empty()?3:(P[u].lab=="C"?0:(P[u].lab=="D"?2:1)));
        b.nm.push_back(P[u].name);
        if (P[u].ch.empty()) { b.leaves.push_back(id); b.sz[id]=1; return id; }
        vector<int> ids; for (int c : P[u].ch) ids.push_back(conv(c));
        if (ids.size()==1) {
            int c=ids[0]; b.lc[id]=b.lc[c]; b.rc[id]=b.rc[c]; b.code[id]=b.code[c]; b.sz[id]=b.sz[c];
            b.nm[id]=b.nm[c];
            if (b.lc[id]>=0){ b.par[b.lc[id]]=id; b.par[b.rc[id]]=id; } else b.leaves.back()=id;
            return id;
        }
        if (ids.size()==2) { b.lc[id]=ids[0]; b.rc[id]=ids[1]; }
        else {
            priority_queue<pair<int,int>, vector<pair<int,int>>, greater<>> pq;
            for (int c : ids) pq.push({b.sz[c], c});
            while (pq.size()>2) {
                auto [s1,a]=pq.top(); pq.pop(); auto [s2,c]=pq.top(); pq.pop();
                int nid=b.n++; b.lc.push_back(a); b.rc.push_back(c); b.par.push_back(-1);
                b.sz.push_back(s1+s2); b.code.push_back(0); b.nm.push_back("");
                b.par[a]=nid; b.par[c]=nid;
                pq.push({s1+s2,nid});
            }
            auto [s1,a]=pq.top(); pq.pop(); auto [s2,c]=pq.top(); pq.pop();
            b.lc[id]=a; b.rc[id]=c;
        }
        b.par[b.lc[id]]=id; b.par[b.rc[id]]=id;
        b.sz[id]=b.sz[b.lc[id]]+b.sz[b.rc[id]];
        return id;
    };
    b.root = conv(root);
    vector<int> st{b.root}; while(!st.empty()){ int t=st.back(); st.pop_back(); b.order.push_back(t); if(b.lc[t]>=0){ st.push_back(b.lc[t]); st.push_back(b.rc[t]); } }
    reverse(b.order.begin(), b.order.end());
    return b;
}

static inline void step1(St& A, const St& U, const St& V, int nu, int nv) {
    auto [Ru,Bu,Gu,RRu,BBu,GGU,RBu,RGu,BGu,coRu,coBu,coGu,Qu] = U;
    auto [Rv,Bv,Gv,RRv,BBv,GGV,RBv,RGv,BGv,coRv,coBv,coGv,Qv] = V;
    A[0]=nu*Ru+nv*Rv; A[1]=nu*Bu+nv*Bv; A[2]=nu*Gu+nv*Gv;
    A[3]=RRu+RRv+nu*nv*Ru*Rv; A[4]=BBu+BBv+nu*nv*Bu*Bv; A[5]=GGU+GGV+nu*nv*Gu*Gv;
    A[6]=RBu+RBv+nu*nv*(Ru*Bv+Bu*Rv); A[7]=RGu+RGv+nu*nv*(Ru*Gv+Gu*Rv); A[8]=BGu+BGv+nu*nv*(Bu*Gv+Gu*Bv);
    A[9]=nu*coRu+nv*coRv + (BBu*nv*Gv + BBv*nu*Gu + GGU*nv*Bv + GGV*nu*Bu + BGu*nv*Rv + BGv*nu*Ru);
    A[10]=nu*coBu+nv*coBv + (RRu*nv*Gv + RRv*nu*Gu + GGU*nv*Rv + GGV*nu*Ru + RGu*nv*Bv + RGv*nu*Bu);
    A[11]=nu*coGu+nv*coGv + (BBu*nv*Rv + BBv*nu*Ru + RRu*nv*Bv + RRv*nu*Bu + RBu*nv*Gv + RBv*nu*Gu);
    A[12]=Qu+Qv + nu*nv*(coRu*Rv+coRv*Ru+coBu*Bv+coBv*Bu+coGu*Gv+coGv*Gu) + (RRu*BGv+RRv*BGu+BBu*RGv+BBv*RGu+GGU*RBv+GGV*RBu);
}
static inline void step2(St& A, const St& U, const St& V, int cu, int cv, int nu, int nv) {
    auto [Ru,Bu,Gu,RRu,BBu,GGU,RBu,RGu,BGu,coRu,coBu,coGu,Qu] = U;
    auto [Rv,Bv,Gv,RRv,BBv,GGV,RBv,RGv,BGv,coRv,coBv,coGv,Qv] = V;
    A[0]=nu*Ru+nv*Rv; A[1]=nu*Bu+nv*Bv; A[2]=nu*Gu+nv*Gv;
    A[3]=cu*RRu+cv*RRv+nu*nv*Ru*Rv; A[4]=cu*BBu+cv*BBv+nu*nv*Bu*Bv; A[5]=cu*GGU+cv*GGV+nu*nv*Gu*Gv;
    A[6]=cu*RBu+cv*RBv+nu*nv*(Ru*Bv+Bu*Rv); A[7]=cu*RGu+cv*RGv+nu*nv*(Ru*Gv+Gu*Rv); A[8]=cu*BGu+cv*BGv+nu*nv*(Bu*Gv+Gu*Bv);
    A[9]=nu*coRu+nv*coRv + cu*nv*(BBu*Gv+GGU*Bv+BGu*Rv) + cv*nu*(BBv*Gu+GGV*Bu+BGv*Ru);
    A[10]=nu*coBu+nv*coBv + cu*nv*(RRu*Gv+GGU*Rv+RGu*Bv) + cv*nu*(RRv*Gu+GGV*Ru+RGv*Bu);
    A[11]=nu*coGu+nv*coGv + cu*nv*(RRu*Bv+BBu*Rv+RBu*Gv) + cv*nu*(RRv*Bu+BBv*Ru+RBv*Gu);
    A[12]=Qu+Qv + nu*nv*(coRu*Rv+coRv*Ru+coBu*Bv+coBv*Bu+coGu*Gv+coGv*Gu) + cu*cv*(RRu*BGv+RRv*BGu+BBu*RGv+BBv*RGu+GGU*RBv+GGV*RBu);
}

struct DP {
    vector<St> V1, V2;
    void init(Bin const& g, vector<int> const& col) {
        V1.assign(g.n, St{}); V2.assign(g.n, St{});
        for (auto& s : V1) s.fill(0);
        for (auto& s : V2) s.fill(0);
        for (int t : g.order) {
            if (g.lc[t] < 0) { int c=col[t]; V1[t].fill(0); V2[t].fill(0); V1[t][c]=1; V2[t][c]=1; continue; }
            int l=g.lc[t], r=g.rc[t];
            int nu=(g.code[l]==2)?0:1, nv=(g.code[r]==2)?0:1, cu=(g.code[l]==0)?1:0, cv=(g.code[r]==0)?1:0;
            step1(V1[t], V1[l], V1[r], nu, nv);
            step2(V2[t], V2[l], V2[r], cu, cv, nu, nv);
        }
    }
    void recolor(Bin const& g, int leaf, int c) {
        V1[leaf].fill(0); V2[leaf].fill(0); V1[leaf][c]=1; V2[leaf][c]=1;
        for (int a=g.par[leaf]; a>=0; a=g.par[a]) {
            int l=g.lc[a], r=g.rc[a];
            int nu=(g.code[l]==2)?0:1, nv=(g.code[r]==2)?0:1, cu=(g.code[l]==0)?1:0, cv=(g.code[r]==0)?1:0;
            step1(V1[a], V1[l], V1[r], nu, nv);
            step2(V2[a], V2[l], V2[r], cu, cv, nu, nv);
        }
    }
};

int main(int argc, char** argv) {
    string spFile, geneFile; int threads=0; bool quiet=false; int limit=0;
    for (int i=1;i<argc;++i) {
        string a=argv[i];
        auto nx=[&]() -> string { if(i+1<argc) return string(argv[++i]); exit(2); };
        if (a=="-t") spFile=nx(); else if (a=="-g") geneFile=nx();
        else if (a=="-j") threads=stoi(nx()); else if (a=="-q") quiet=true;
        else if (a=="--limit") limit=stoi(nx());
        else { cerr<<"unknown "<<a<<"\n"; return 2; }
    }
    vector<string> spL, gL;
    { ifstream f(spFile); string l; while(getline(f,l)){ bool nz=false; for(char c:l) if(!isspace((unsigned char)c)){nz=true;break;} if(nz) spL.push_back(l);} }
    { ifstream f(geneFile); string l; while(getline(f,l)){ bool nz=false; for(char c:l) if(!isspace((unsigned char)c)){nz=true;break;} if(nz) gL.push_back(l);} }
    if (spL.empty() || gL.empty()) { cerr<<"empty input\n"; return 1; }
    if (limit>0 && (int)gL.size()>limit) gL.resize(limit);

    int spRoot; vector<Raw> spP = parseRaw(spL[0], spRoot);
    // 物种树：根三叉 -> 合并后两个孩子
    if (spP[spRoot].ch.size()==3) {
        int synth = spP.size(); spP.push_back({}); spP[synth].ch = { spP[spRoot].ch[1], spP[spRoot].ch[2] };
        spP[spRoot].ch.resize(1); spP[spRoot].ch.push_back(synth);
    }
    // 二叉检查
    for (auto& nd : spP) if (!nd.ch.empty() && nd.ch.size()!=2) { cerr<<"species tree not binary\n"; return 1; }

    unordered_map<string,int> name2id; vector<string> id2name;
    for (auto& nd : spP) if (nd.ch.empty() && !name2id.count(nd.name)) { name2id[nd.name]=id2name.size(); id2name.push_back(nd.name); }
    int nTaxa = id2name.size();
    int W = (nTaxa+63)/64;
    vector<vector<u64>> bits(spP.size(), vector<u64>(W,0));
    vector<int> ord; { vector<pair<int,int>> st{{spRoot,0}}; while(!st.empty()){ auto [u,ph]=st.back(); st.pop_back();
        if(!ph){ st.push_back({u,1}); for(auto it=spP[u].ch.rbegin(); it!=spP[u].ch.rend(); ++it) st.push_back({*it,0}); }
        else ord.push_back(u); } }
    for (int u : ord) { if (spP[u].ch.empty()) { int id=name2id[spP[u].name]; bits[u][id>>6]|=1ULL<<(id&63); }
        else for (int c : spP[u].ch) for (int w=0;w<W;++w) bits[u][w]|=bits[c][w]; }
    vector<int> spNodes; for (int u : ord) if (!spP[u].ch.empty()) spNodes.push_back(u);

    int T=(int)gL.size(), P=threads>0?threads:omp_get_max_threads();
    vector<ll> perGene(T,0); vector<char> skipped(T,0);
    auto t0=chrono::steady_clock::now();
    #pragma omp parallel for schedule(dynamic) num_threads(P)
    for (int gi=0; gi<T; ++gi) {
        int gRoot; vector<Raw> gP = parseRaw(gL[gi], gRoot);
        int nLeaf=0; for (auto& nd : gP) if (nd.ch.empty()) ++nLeaf;
        if (nLeaf<4) { skipped[gi]=1; continue; }
        Bin g = toBinary(gP, gRoot);
        vector<int> tid(g.n, -1);
        for (int lf : g.leaves) { auto it=name2id.find(g.nm[lf]); if (it!=name2id.end()) tid[lf]=it->second; }
        vector<int> col(g.n, 2);
        DP dp; dp.init(g, col);
        ll sum1=0, sum2=0;
        // 逐物种节点：从零染色（验证语义用；增量版本见 recolor）
        for (int u : spNodes) {
            int LC=spP[u].ch[0], RC=spP[u].ch[1];
            // 三划分：色0 = RC 的叶子，色1 = LC 的叶子，色2 = 其余
            u64 const* A = bits[RC].data();
            vector<u64> Bbits = bits[LC];
            for (int lf : g.leaves) {
                int t=tid[lf]; int c=2;
                if (t>=0) {
                    if ((A[t>>6]>>(t&63))&1ULL) c=0;
                    else if ((Bbits[t>>6]>>(t&63))&1ULL) c=1;
                }
                if (col[lf]!=c) { col[lf]=c; dp.recolor(g, lf, c); }
            }
            if (getenv("DBG")) {
                fprintf(stderr, "[dbg] species node %d: LC=%d RC=%d  colors:", u, LC, RC);
                for (int lf : g.leaves) fprintf(stderr, " %s=%d", g.nm[lf].c_str(), col[lf]);
                fprintf(stderr, "  V1root=%lld V2root=%lld\n", dp.V1[g.root][12], dp.V2[g.root][12]);
            }
            if (getenv("DBG2")) { fprintf(stderr, "[mine] node %d (LC=%d RC=%d):", u, LC, RC);
                for (int lf : g.leaves) fprintf(stderr, " %s(t=%d)=%d", g.nm[lf].c_str(), tid[lf], col[lf]);
                fprintf(stderr, "  Q1=%lld Q2=%lld\n", dp.V1[g.root][12], dp.V2[g.root][12]); }
            sum1 += dp.V1[g.root][12];
            sum2 += dp.V2[g.root][12];
        }
        perGene[gi] = sum1/2 - sum2/2;

    }
    double secs = chrono::duration<double>(chrono::steady_clock::now()-t0).count();
    __int128 total=0; int nSkip=0; for (int i=0;i<T;++i){ total+=perGene[i]; if(skipped[i]) ++nSkip; }
    auto i128str=[]( __int128 v){ if(!v) return string("0"); string s; bool neg=v<0; unsigned __int128 u=neg?(unsigned __int128)(-v):(unsigned __int128)v; while(u){s.push_back('0'+u%10); u/=10;} if(neg) s.push_back('-'); reverse(s.begin(),s.end()); return s; };
    cout << "{\"score\":" << i128str(total) << ",\"n_taxa\":" << nTaxa << ",\"n_genes\":" << T
         << ",\"n_scored_genes\":" << (T-nSkip) << ",\"n_skipped_genes\":" << nSkip
         << ",\"runtime_sec\":" << fixed << setprecision(3) << secs << ",\"per_gene\":[";
    for (int i=0;i<T;++i){ if(i) cout<<','; cout<<perGene[i]; }
    cout << "]}\n";
    return 0;
}
