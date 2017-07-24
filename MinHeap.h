#pragma once

#include<vector>
#include<algorithm>
#include<limits>

template <typename T>
class MinHeap
{
    using sz_t = unsigned long long;
    using vl_t = T;
    const sz_t INVLID = std::numeric_limits<sz_t>::max();

    sz_t max_size;
    std::vector<std::pair<sz_t,T>> tree;
    std::vector<sz_t> id_to_pos;
    void swap_node(int a,int b)
    {
        std::swap( id_to_pos[tree[a].first] , id_to_pos[tree[b].first] );
        std::swap( tree[a] , tree[b] );
    }

    sz_t father(sz_t id)
    {
        return (id-1)/2;
    }

    sz_t lchild(sz_t id)
    {
        return id*2+1;
    }

    sz_t rchild(sz_t id)
    {
        return id*2+2;
    }

    bool idok(sz_t id)
    {
        return id < tree.size();
    }

    void up(sz_t id)
    {
        while( id && tree[father(id)].second > tree[id].second )
        {
            swap_node(id,father(id));
            id = father(id);
        }
    }

    void down(sz_t id)
    {
        sz_t target;
        while( idok( target = lchild(id) ) )
        {
            if( idok( rchild(id) ) && tree[target].second > tree[rchild(id)].second )
                target = rchild(id);
            if( tree[id].second > tree[target].second )
            {
                swap_node(id,target);
                id = target;
            }
            else break;
        }
    }
public:

    MinHeap(sz_t max_sz)
    {
        max_size = max_sz;
        id_to_pos.resize(max_size,INVLID);
    }

    void push(sz_t index,vl_t val)
    {
        //if index exist, use update
        if( id_to_pos[index]!=INVLID )
        {
            update(index,val);
            return ;
        }
        id_to_pos[index] = tree.size();
        tree.emplace_back(index,val);
        up(id_to_pos[index]);
    }

    void update(sz_t index,vl_t val)
    {
        sz_t tid = id_to_pos[index];
        bool isless = val < tree[tid].second ? true : false;
        tree[tid].second = val;

        if( isless ) up(tid);
        else  down(tid);
    }

    std::pair<sz_t,vl_t> top()
    {
        return tree[0];
    }

    void pop()
    {
        id_to_pos[ tree[0].first ] = INVLID;
        tree[0] = tree.back();

        tree.pop_back();
        if( !tree.empty() )
            down(0);
    }

    void clear()
    {
        for(auto &v:id_to_pos) v = INVLID;
        tree.clear();
    }

    bool empty()
    {
        return tree.empty();
    }

    u64 size()
    {
        return tree.size();
    }
};
