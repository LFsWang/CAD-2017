#pragma once

#include<vector>
#include<algorithm>
/*
    DisjoinSet use 0-base
*/
class DisjoinSet
{
public:
    DisjoinSet(int=0);
    void init(int);
    int find(int);
    bool same(int,int);
    void U(int,int);
    int size(int);
private:
    int sz;
    std::vector<int> dis;
    std::vector<int> sum;
};

DisjoinSet::DisjoinSet(int s)
{
    init(s);
}

void DisjoinSet::init(int s)
{
    sz=s;
    dis.resize(sz);
    sum.resize(sz);
    for(int i=0;i<sz;++i)
    {
        dis[i]=i;
        sum[i]=1;
    }
}

int DisjoinSet::find(int p)
{
    if( dis[p]==p )
        return p;
    return dis[p]=find(dis[p]);
}

bool DisjoinSet::same(int a,int b)
{
    return find(a)==find(b);
}

void DisjoinSet::U(int a,int b)
{
    if( !same(a,b) )
    {
        if( sum[dis[a]] < sum[dis[b]] )
            std::swap(a,b);
        sum[dis[a]] += sum[dis[b]];
        dis[dis[b]] = dis[a];
    }
}

int DisjoinSet::size(int p)
{
    return sum[find(p)];
}

