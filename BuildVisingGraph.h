#pragma once
#include<map>
#include<set>
#include<tuple>
#include<algorithm>

struct Edge{
	u32 ori_u;
	u32 ori_v;
	u8 type;
	//建樹只需要u,v,cost
	u32 u;
	u32 v;
	u64 cost;
};

struct point3D{
	s64 x;
	s64 y;
	u32 layer;
	bool operator<(const point3D &b)const{
		if(x!=b.x)return x<b.x;
		if(y!=b.y)return y<b.y;
		return layer<b.layer;
	}
};

struct VisingGraph{
	u32 N;
	std::vector<std::vector<u32>> G;
	std::vector<Edge> edge;
	std::vector<bool> is_pinv;
	//建樹只需要N,G,edge,is_pinv
	
	std::vector<std::vector<u32>> ori_G;
	std::vector<std::pair<point3D,u8>> V_set;
	std::vector<u32> shrink_from;
	std::vector<s64> Px;
	std::vector<s64> Py;
	
	void build(const DataSet &data);
};