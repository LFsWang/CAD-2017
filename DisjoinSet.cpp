#include"DisjoinSet.h"

DisjoinSet::DisjoinSet(size_t s)
{
    init(s);
}

void DisjoinSet::init(size_t s)
{
    sz=s;
    dis.resize(sz);
    sum.resize(sz);
    for(size_t i=0;i<sz;++i)
    {
        dis[i]=i;
        sum[i]=1;
    }
}

size_t DisjoinSet::find(size_t p)
{
    if( dis[p]==p )
        return p;
    return dis[p]=find(dis[p]);
}

bool DisjoinSet::same(size_t a,size_t b)
{
    return find(a)==find(b);
}

void DisjoinSet::U(size_t a,size_t b)
{
    if( !same(a,b) )
    {
        if( sum[dis[a]] < sum[dis[b]] )
            std::swap(a,b);
        sum[dis[a]] += sum[dis[b]];
        dis[dis[b]] = dis[a];
    }
}

size_t DisjoinSet::size(size_t p)
{
    return sum[find(p)];
}