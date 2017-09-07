#include "DataLoader.h"
#include "BuildVisingGraph.h"
#include "DisjoinSet.h"
#include "Queue.h"

#include<ctime>
#include<functional>
#include<vector>
#include<cassert>
#include<iostream>
#include<cstdio>
#include<tuple>
#include<unordered_set>
#include<stack>

using std::cout;
using std::endl;

inline void showclock(const char *str=nullptr)
{
#ifdef _WIN32
    long long CL_PER_SEC = CLOCKS_PER_SEC;//1000;
#else
    long long CL_PER_SEC = CLOCKS_PER_SEC;//1000000;//test on centos
#endif
    static long long last = 0;
    auto show_time = [&](long long time){
        long long ms = time%CL_PER_SEC; time/=CL_PER_SEC;
        long long sec = time%60;  time/=60;
        long long min = time%60;  time/=60;
        printf("%2llu:%02llu:%02llu %06llu",time,min,sec,ms);
    };
    
    long long now = std::clock();
    if(str)printf("%s ,",str);
    printf("Time:");show_time(now);printf("\t(");
    long long diff = now - last;
    show_time(diff);printf(")\n");
    last = now;
}

std::vector<std::size_t> select_edge(const VisingGraph &G,DisjoinSet &ds)
{
    using sz_t = std::size_t;
    const sz_t INVLID = std::numeric_limits<sz_t>::max();
    const u64 INF = 0x3fffffffffffffffLL;
    static_assert( INF <= std::numeric_limits<decltype(INF)>::max()/2 ,"Invlid inf!");

    const sz_t N = G.G.size();
    auto &V = G.G;

    std::vector<u64>  dist;
    std::vector<sz_t> prev_eid;
    std::vector<sz_t> index;

    //SPFA
    std::vector<bool> inqueue(N,false);
    Queue<sz_t> qu;
    dist.reserve(N);
    prev_eid.reserve(N);
    index.reserve(N);
    for(sz_t i=0;i<N;++i)
    {
        dist.emplace_back(INF);
        prev_eid.emplace_back(INVLID);
        index.emplace_back(INVLID);
        if( G.is_pinv[i] )
        {
            dist[i] = 0;
            index[i]=i;
            qu.push_back(i);
            inqueue[i]=true;
        }
    }
    showclock(" :INIT FIND EDGE");
    sz_t v;
    while( !qu.empty() )
    {
        v = qu.front();
        qu.pop_front();
        inqueue[v]=false;
        for(sz_t eid:V[v])
        {
            sz_t e = G.edge[eid].v;
            u64 cost = G.edge[eid].cost;

            if( dist[e] > dist[v]+cost )
            {
                dist[e] = dist[v]+cost;
                prev_eid[e] = eid;
                index[e] = index[v];
                if( !inqueue[e] )
                {
                    if( !qu.empty() && dist[e] <= dist[qu.front()] )
                        qu.push_front(e);
                    else
                        qu.push_back(e);
                    inqueue[e]=true;
                }
            }
        }
    }
    showclock(" :SPFA");

    std::vector< std::tuple<u64,sz_t> > CrossEdge;
    {
        sz_t eid=0;
        for(const Edge &E:G.edge)
        {
            if( index[E.u]!=index[E.v] && E.u > E.v )
            {
                CrossEdge.emplace_back( dist[E.u]+dist[E.v]+E.cost , eid );
            }
            eid++;
        }
    }
    std::sort(CrossEdge.begin(),CrossEdge.end());
    showclock(" :Find CrossEdge");

    std::vector<sz_t> &SelectKEdge = qu.shift();
    SelectKEdge.clear();
    ds.init(N);
    for(const auto &TUS:CrossEdge)
    {
        sz_t eid;
        std::tie(std::ignore,eid) = TUS;
        const auto &E = G.edge[eid];
        if( !ds.same(index[E.u],index[E.v]) )
        {
            SelectKEdge.emplace_back(eid);
            ds.U(index[E.u],index[E.v]);
        }
    }
    showclock(" :Kruskal 1");

    std::vector<bool> &used = inqueue;
    used.resize(G.edge.size());
    std::fill(used.begin(),used.end(),false);
    
    auto &FinalEdge = CrossEdge;
    FinalEdge.clear();
    for(sz_t eid:SelectKEdge)
    {
        FinalEdge.emplace_back(G.edge[eid].cost,eid);
        for(int t=0;t<2;++t)
        {
            sz_t v = (t&1)?G.edge[eid].u:G.edge[eid].v;
            while( v!=index[v] && !used[prev_eid[v]] )
            {
                used[prev_eid[v]] = true;
                FinalEdge.emplace_back(G.edge[prev_eid[v]].cost,prev_eid[v]);
                v = G.edge[prev_eid[v]].u;
            }
        }
    }
    showclock(" :Recover Edge");

    auto &deg = dist;
    std::fill(deg.begin(),deg.end(),0);
    std::unordered_set<int> used_eid;
    
    std::sort(FinalEdge.begin(),FinalEdge.end());
    ds.init(N);
    SelectKEdge.clear();
    for(const auto &TUS:CrossEdge)
    {
        sz_t eid;
        std::tie(std::ignore,eid) = TUS;
        const auto &E = G.edge[eid];
        if( !ds.same(E.u,E.v) )
        {
            used_eid.insert(eid);
            deg[E.u]++;
            deg[E.v]++;
            SelectKEdge.emplace_back(eid);
            ds.U(E.u,E.v);
        }
    }
    showclock(" :Kruskal 2");

    std::stack<sz_t> stack;
    for(sz_t i=0;i<N;++i)
    {
        if( !G.is_pinv[i] && deg[i]==1 )
        {
            deg[i]--;
            stack.emplace(i);
        }
    }

    std::cout<<"Before reduce E="<<used_eid.size()<<",stack hold:"<<stack.size()<<std::endl;
    while( !stack.empty() )
    {
        sz_t v = stack.top();
        stack.pop();

        //find the edge to delete
        for(sz_t eid:V[v])
        {
            if( used_eid.find(eid) == used_eid.end() )
                continue;
            used_eid.erase(eid);

            auto e = G.edge[eid].v;
            deg[ e ]--;
            if( deg[ e ]==0 && !G.is_pinv[e] )
                stack.emplace(e);
        }
    }
    SelectKEdge.clear();
    std::copy(used_eid.begin(),used_eid.end(),std::back_inserter(SelectKEdge));
    std::cout<<"After reduce E="<<SelectKEdge.size()<<std::endl;
    showclock(" :Reduce Edge");

    return SelectKEdge;
}

int main(int argc,char *argv[])
{
    //show compile info
    std::cout<<"JINKELA NET_OPEN_FINDER"<<std::endl;
    std::cout<<"Compile time:"<<__DATE__<<' '<<__TIME__<<std::endl;
#ifdef __GNUC__
    std::cout<<"G++ version: "<<__VERSION__<<endl;
#endif
    std::cout<<"======================================"<<endl;
    DataSet d;
    std::ifstream fin;
    std::ofstream fout;
    if( argc>1 )
        fin.open(argv[1]);
    else
        fin.open("a.in");
    
    if( argc>2 )
        fout.open(argv[2]);
    else
        fout.open("ans.out");
    if( !fin.is_open() )
    {
        std::cout<<"Open Input File fail!"<<std::endl;
        exit(-1);
    }
    d.load( fin );
    showclock("Load File");

	d.set_spacing_on_Obstacles();
    showclock("set_spacing_on_Obstacles");

	VisingGraph v;
	v.build(d,1);
    showclock("VisingGraph build");

	std::vector<std::size_t> res=select_edge(v,v.DST);
    showclock("select_edge");
	
    v.print_select_edges(res,fout);
	
    showclock("DONE!!!");
}
