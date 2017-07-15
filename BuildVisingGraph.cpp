#include"DataLoader.h"
#include"BuildVisingGraph.h"
#include"SwapLineP1.h"

inline void set_dis(std::vector<s64> &dis)
{
	std::sort(dis.begin(),dis.end());
	dis.resize(std::unique(dis.begin(),dis.end())-dis.begin());
}
inline u32 get_dis(const std::vector<s64> &dis,s64 val)
{
	return std::lower_bound(dis.begin(),dis.end(),val)-dis.begin();
}
inline u32 get_upper_dis(const std::vector<s64> &dis,s64 val)
{
	return std::upper_bound(dis.begin(),dis.end(),val)-dis.begin();
}

inline void get_ori_P(std::vector<s64> &Px,std::vector<s64> &Py,const DataSet &data)
{
	for(s32 i=1;i<=data.metal_layers;++i)
	{
		for(auto o:data.Obstacles[i])
		{
			Px.emplace_back(o.first.x);
			Px.emplace_back(o.second.x);
			Py.emplace_back(o.first.y);
			Py.emplace_back(o.second.y);
		}
		for(auto s:data.RoutedShape[i])
		{
			Px.emplace_back(s.first.x);
			Px.emplace_back(s.second.x);
			Py.emplace_back(s.first.y);
			Py.emplace_back(s.second.y);
		}
	}
	set_dis(Px);
	set_dis(Py);
}

inline void get_Line(std::vector<std::tuple<s64,s64,s64>> &line,std::vector<statementP1> &state,swape_line_P1 &SL,u32 l,u32 r)
{
	sort(state.begin(),state.end());
	SL.init(r-l+1);
	
	for(const auto &st:state)
	{
		if(st.inout>0)
		{
			SL.find_seg(st.b1,st.b2,l,r,1);
			SL.insert(st.b1,st.b2,1,l,r,1);
		}
		else
		{
			SL.insert(st.b1,st.b2,-1,l,r,1);
			SL.find_seg(st.b1,st.b2,l,r,1);
		}
		SL.set_seg();
		for(const auto &seg:SL.segments)
		{
			line.emplace_back(st.a,seg.first,seg.second);
		}
		SL.segments.clear();
	}
}

inline void get_xyLine(s32 lay,std::vector<std::tuple<s64,s64,s64>> &xLine,std::vector<std::tuple<s64,s64,s64>> &yLine,const DataSet &data,const std::vector<s64> &Px,const std::vector<s64> &Py)
{
	std::vector<statementP1> state;
	swape_line_P1 SL;
	
	for(auto o:data.Obstacles[lay])
	{
		u32 x1=get_dis(Px,o.first.x);
		u32 x2=get_dis(Px,o.second.x);
		u32 y1=get_dis(Py,o.first.y);
		u32 y2=get_dis(Py,o.second.y);
		
		state.emplace_back(x1,y1,y2, 1);
		state.emplace_back(x2,y1,y2,-1);
	}
	for(auto o:data.RoutedShape[lay])
	{
		u32 x1=get_dis(Px,o.first.x);
		u32 x2=get_dis(Px,o.second.x);
		u32 y1=get_dis(Py,o.first.y);
		u32 y2=get_dis(Py,o.second.y);
		
		state.emplace_back(x1,y1,y2, 1);
		state.emplace_back(x2,y1,y2,-1);
	}
	
	get_Line(xLine,state,SL,0,Py.size());
	
	
	state.clear();
	for(auto o:data.Obstacles[lay])
	{
		u32 x1=get_dis(Px,o.first.x);
		u32 x2=get_dis(Px,o.second.x);
		u32 y1=get_dis(Py,o.first.y);
		u32 y2=get_dis(Py,o.second.y);
		
		state.emplace_back(y1,x1,x2, 1);
		state.emplace_back(y2,x1,x2,-1);
	}
	for(auto o:data.RoutedShape[lay])
	{
		u32 x1=get_dis(Px,o.first.x);
		u32 x2=get_dis(Px,o.second.x);
		u32 y1=get_dis(Py,o.first.y);
		u32 y2=get_dis(Py,o.second.y);
		
		state.emplace_back(y1,x1,x2, 1);
		state.emplace_back(y2,x1,x2,-1);
	}
	
	get_Line(yLine,state,SL,0,Px.size());
	
	for(auto i:xLine){
		std::get<0>(i)=Px[std::get<0>(i)];
		std::get<1>(i)=Py[std::get<1>(i)];
		std::get<2>(i)=Py[std::get<2>(i)];
	}
	for(auto i:yLine){
		std::get<0>(i)=Py[std::get<0>(i)];
		std::get<1>(i)=Px[std::get<1>(i)];
		std::get<2>(i)=Px[std::get<2>(i)];
	}
}

inline void reGet_ori_P(std::vector<s64> &Px,std::vector<s64> &Py,const DataSet &data,std::vector<std::tuple<s64,s64,s64>> *xLine,std::vector<std::tuple<s64,s64,s64>> *yLine)
{
	Px.clear();
	Py.clear();
	
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		for(const auto &i:xLine[lay])
		{
			Px.emplace_back(std::get<0>(i));
			Py.emplace_back(std::get<1>(i));
			Py.emplace_back(std::get<2>(i));
		}
		
		for(const auto &i:yLine[lay])
		{
			Py.emplace_back(std::get<0>(i));
			Px.emplace_back(std::get<1>(i));
			Px.emplace_back(std::get<2>(i));
		}
		for(const auto &i:data.RoutedVia[lay])
		{
			Px.emplace_back(i.x);
			Py.emplace_back(i.y);
		}
	}
	Px.emplace_back(data.boundary.first.x);
	Px.emplace_back(data.boundary.second.x);
	Py.emplace_back(data.boundary.first.y);
	Py.emplace_back(data.boundary.second.y);
	
	set_dis(Px);
	set_dis(Py);
}

void set_2d_VG_point(s32 l,s32 r)
{
	if(l==r)return;
	s32 mid=(l+r)/2;
	
}

void VisingGraph::build(const DataSet &data)
{
	static const int LIMIT_LAYER = DataSet::LIMIT_LAYER;
	
	std::vector<std::tuple<s64,s64,s64>> xLine[LIMIT_LAYER];
	std::vector<std::tuple<s64,s64,s64>> yLine[LIMIT_LAYER];
	
	get_ori_P(Px,Py,data);
	
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		get_xyLine(lay,xLine[lay],yLine[lay],data,Px,Py);
	}
	
	reGet_ori_P(Px,Py,data,xLine,yLine);
	
	std::vector<std::pair<u32,u32>> P1[LIMIT_LAYER];
	
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		for(const auto &line:xLine[lay])
		{
			u32 x=get_dis(Px,std::get<0>(line));
			u32 l=get_dis(Py,std::get<1>(line));
			u32 r=get_dis(Py,std::get<2>(line));
			for(;l<=r;++l)P1[lay].emplace_back(x,l);
		}
		for(const auto &line:yLine[lay])
		{
			u32 y=get_dis(Py,std::get<0>(line));
			u32 l=get_dis(Px,std::get<1>(line));
			u32 r=get_dis(Px,std::get<2>(line));
			for(;l<=r;++l)P1[lay].emplace_back(l,y);
		}
	}
	
	
}