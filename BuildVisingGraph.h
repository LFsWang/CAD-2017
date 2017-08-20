#pragma once
#include<map>
#include<set>
#include<tuple>
#include<algorithm>

struct Edge{
	size_t ori_u;
	size_t ori_v;
	u8 type;
	//建樹只需要u,v,cost
	size_t u;
	size_t v;
	u64 cost;
	Edge(size_t ori_u,size_t ori_v,u8 type):ori_u(ori_u),ori_v(ori_v),type(type){}
};

struct point3D{
	u32 x;
	u32 y;
	u32 layer;
	point3D(){}
	point3D(u32 x,u32 y,u32 z):x(x),y(y),layer(z){}
	bool operator<(const point3D &b)const{
		if(layer!=b.layer)return layer<b.layer;
		if(x!=b.x)return x<b.x;
		return y<b.y;
	}
	bool operator==(const point3D &b)const{
		return x==b.x&&y==b.y&&layer==b.layer;
	}
};

struct VisingGraph{
	size_t N;
	std::vector<std::vector<size_t>> G;
	std::vector<Edge> edge;
	std::vector<bool> is_pinv;
	//建樹只需要N,G,edge,is_pinv
	
	std::vector<std::vector<u32>> ori_G;
	std::vector<point3D> V_set;
	std::vector<size_t> shrink_from;
	std::vector<s64> Px;
	std::vector<s64> Py;
	
	void build(const DataSet &data,bool is_not_connect);
};