#pragma once
#include<map>
#include<set>
#include<tuple>
#include<iostream>
#include<algorithm>
#include<unordered_map>
#include"DisjoinSet.h"
#include"SwapLineStatement.h"

struct Edge{
	size_t ori_u;
	size_t ori_v;
	u8 type;
	//建樹只需要u,v,cost
	size_t u;
	size_t v;
	u32 cost;
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
	
	std::vector<point3D> V_set;
	std::vector<size_t> shrink_from;
	std::vector<s64> Px;
	std::vector<s64> Py;
	std::vector<Statemant_2D_VG> xLine[10+1];
	std::vector<Statemant_2D_VG> yLine[10+1];
	
	DisjoinSet DST;
	
	void build(const DataSet &data,bool is_not_connect);
	void build_beta(const DataSet &data,bool is_not_connect);
	void print_select_edges(const std::vector<std::size_t> &res,std::ofstream &fout);
};

struct pairhash {
	template <typename T, typename U>
	std::size_t operator()(const std::pair<T,U> &x) const
	{
		size_t seed = std::hash<T>()(x.first);
		return std::hash<U>()(x.second) + 0x9e3779b9 + (seed<<6) + (seed>>2);
		//return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
	}
};