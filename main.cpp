#include "DataLoader.h"
#include "BuildVisingGraph.h"
#include "DisjoinSet.h"
#include "MinHeap.h"

#include<ctime>
#include<functional>
#include<vector>
#include<cassert>
#include<iostream>
#include<cstdio>

using std::cout;
using std::endl;
std::vector<u64> find_steiner_point_id(const VisingGraph &G)
{
    const u64 INVLID = std::numeric_limits<u64>::max();
    const s64 INF = 0x3fffffffffffffffLL;
    u64 N = G.G.size();
    auto &V = G.G;

    std::vector<s64> dist(N,INF);
    std::vector<u64> stack;
    std::vector<u64> steiner_point_id;
    DisjoinSet ds(N);

    using min_heap = MinHeap<s64>;
    min_heap pq(N);
//cout<<"N="<<N<<'\n';
    for(u64 i=0;i<N;++i)
    {
        if( G.is_pinv[i] )
        {
            pq.push(i,0);
            dist[i] = 0;
            //cout<<"ping:"<<i<<endl;
            cout<<pq.size()<<endl;
        }
    }

    for(u64 T=0;T<N-1;++T)
    {
        u64 traget = INVLID;
        while( !pq.empty() )
        {
            u64 v,e;
            s64 cost;
            std::tie(v,cost) = pq.top();
            pq.pop();//cout<<"view"<<v<<' '<<cost<<endl;
            //promise no 0-w edge between pin point
            if( dist[v]!=0 && G.is_pinv[v] )
            {
                traget = v;
                break;
            }
            for( u64 eid:V[v] )
            {
                e = G.edge[eid].v;
                cost = G.edge[eid].cost;
                if( ds.same(v,e) )//inside merged conpoment
                    cost = 0;
                if( dist[e] > dist[v]+cost )
                {
                    dist[e] = cost;
                    pq.push(e,cost);
                }
            }
        }//cout<<"N="<<traget<<'\n';
        assert(traget!=INVLID);
        //merge sp path to source
        stack.clear();
        stack.emplace_back(traget);
        while( !stack.empty() )
        {
            u64 v = stack.back();
            stack.pop_back();
            if( !G.is_pinv[v] )// paper 3-3
                steiner_point_id.push_back(v);

            if( dist[v]==0 )//reach face of conpoment
                continue;

            for( u64 eid:V[v] )
            {
                u64 e = G.edge[eid].v;
                s64 cost = G.edge[eid].cost;
                if( dist[v] == dist[e]+cost && !ds.same(v,e) )
                {
                    ds.U(v,e);
                    stack.emplace_back(e);
                }
            }
        }
        //reset pq
        pq.clear();
        for(u64 i=0;i<N;++i)
        {
            if( G.is_pinv[i] || ds.size(i)!=1 )
            {
                dist[i] = 0;
                pq.push(i,0);
            }
            else
            {
                dist[i] = INF;
            }
        }
    }
    return steiner_point_id;
}

inline void showclock(const char *str=nullptr)
{
    u64 time = std::clock();
    u64 ms = time%1000; time/=1000;
    u64 sec = time%60;  time/=60;
    u64 min = time%60;  time/=60;
    if(str)printf("%s ,",str);
    printf("Time: %2llu:%02llu:%02llu %03llu\n",time,min,sec,ms);
}

int main(int argc,char *argv[])
{
    DataSet d;
    std::ifstream fin;
    if( argc>1 )
        fin.open(argv[1]);
    else
        fin.open("a.in");
    
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
	v.build(d);
    showclock("VisingGraph build");

	std::vector<u64> res=find_steiner_point_id(v);
    showclock("find_steiner_point_id");
}
