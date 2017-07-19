#pragma once

#include<vector>
#include<algorithm>
/*
    DisjoinSet use 0-base
*/
class DisjoinSet
{
public:
    DisjoinSet(size_t=0);
    void init(size_t);
    size_t find(size_t);
    bool same(size_t,size_t);
    void U(size_t,size_t);
    size_t size(size_t);
private:
    size_t sz;
    std::vector<size_t> dis;
    std::vector<size_t> sum;
};

