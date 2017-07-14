#include "DataLoader.h"
#include "BuildVisingGraph.h"
#include "DisjoinSet.h"

#include<functional>
#include<vector>
#include<queue>
#include<cassert>
std::vector<int> find_steiner_point_id(const VisingGraph &G)
{
    const s64 INF = 0x3fffffffffffffffLL;
    s64 N = G.G.size();
    auto &V = G.G;

    std::vector<s64> dist(N,INF);
    std::vector<s32> stack;
    std::vector<s32> steiner_point_id;
    DisjoinSet ds(N);

    using node = std::pair<s64,s32>;
    using min_heap = std::priority_queue<node,std::vector<node>,std::greater<node>>;
    min_heap pq;

    for(s32 i=0;i<N;++i)
    {
        if( G.is_pinv[i] )
        {
            pq.emplace(0,i);
            dist[i] = 0;
        }
    }

    for(s32 T=0;T<N-1;++T)
    {
        s32 traget = -1;
        while( !pq.empty() )
        {
            s32 v,e;
            s64 cost;
            std::tie(std::ignore,v) = pq.top();
            pq.pop();
            //promise no 0-w edge between pin point
            if( dist[v]!=0 && G.is_pinv[v] )
            {
                traget = v;
                break;
            }
            for( s32 eid:V[v] )
            {
                e = G.edge[eid].v;
                cost = G.edge[eid].cost;
                if( ds.same(v,e) )//inside merged conpoment
                    cost = 0;
                if( dist[e] > dist[v]+cost )
                {
                    dist[e] = cost;
                    pq.emplace(cost,e);
                }
            }
        }
        assert(traget!=-1);
        //merge sp path to source
        stack.clear();
        stack.emplace_back(traget);
        while( !stack.empty() )
        {
            s32 v = stack.back();
            stack.pop_back();
            if( !G.is_pinv[v] )// paper 3-3
                steiner_point_id.push_back(v);

            if( dist[v]==0 )//reach face of conpoment
                continue;

            for( s32 eid:V[v] )
            {
                s32 e = G.edge[eid].v;
                s64 cost = G.edge[eid].cost;
                if( dist[v] == dist[e]+cost && !ds.same(v,e) )
                {
                    ds.U(v,e);
                    stack.emplace_back(e);
                }
            }
        }
        //reset pq
        pq = min_heap();//hack for clear
        for(s32 i=0;i<N;++i)
        {
            if( G.is_pinv[i] || ds.size(i)!=1 )
            {
                dist[i] = 0;
                pq.emplace(0,i);
            }
            else
            {
                dist[i] = INF;
            }
        }
    }
    return steiner_point_id;
}

int main()
{
    DisjoinSet ds(87);
    DataSet d;
    std::ifstream fin("a.in");
    d.load( fin );
	d.set_spacing_on_Obstacles();
}
