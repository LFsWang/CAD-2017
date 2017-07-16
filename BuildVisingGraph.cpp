#include"DataLoader.h"
#include"BuildVisingGraph.h"
#include"SwapLineP1.h"
#include"SwapLineStatement.h"
#include"BinaryIndexTree.h"

template<typename T>
inline void set_dis(T &dis)
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
	
	for(auto &i:xLine){
		std::get<0>(i)=Px[std::get<0>(i)];
		std::get<1>(i)=Py[std::get<1>(i)];
		std::get<2>(i)=Py[std::get<2>(i)];
	}
	for(auto &i:yLine){
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

inline void set_3d_VG_point(std::vector<std::pair<u32,u32>> &S,s32 la1,s32 la2,const DataSet &data,std::vector<std::pair<u32,u32>> *P1,std::vector<s64> &Px,std::vector<s64> &Py)
{
	std::vector<Statemant_2D_VG> state;
	
	for(const auto &i:P1[la1]) S.emplace_back(i);
	
	for(const auto &p:S)
	{
		state.emplace_back(2,p.first,p.second);
	}
	
	for(const auto &o:data.Obstacles[la2])
	{
		s32 x1=get_dis(Px,o.first.x);
		s32 x2=s32(get_upper_dis(Px,o.second.x))-1;
		s32 y1=get_dis(Py,o.first.y);
		s32 y2=get_upper_dis(Py,o.second.y);
		
		if(x1<x2&&y1<y2)
		{
			state.emplace_back(3,x1,y1,y2);
			state.emplace_back(1,x2,y1,y2);
		}
	}
	
	sort(state.begin(),state.end());
	
	BIT ST;
	ST.init(Py.size()+2);
	
	std::vector<std::pair<u32,u32>> nS;
	
	for(const auto &st:state)
	{
		switch(st.type)
		{
			case 1:{
				ST.add(st.b1+1,-1);
				ST.add(st.b2+1,1);
				break;
			}
			case 2:{
				if(ST.get_sum(st.b1+1)==0)
				{
					nS.emplace_back(st.a,st.b1);
				}
				break;
			}
			case 3:{
				ST.add(st.b1+1,1);
				ST.add(st.b2+1,-1);
				break;
			}
			default:
				cout<<"ERROR!\n";
		}
	}
	
	S.swap(nS);
}

void recursive_set_3d_VG_point(s32 l,s32 r,const DataSet &data,std::vector<std::pair<u32,u32>> *P1,std::vector<s64> &Px,std::vector<s64> &Py)
{
	if(l>=r)return;
	s32 mid=(l+r)/2;
	
	std::vector<std::pair<u32,u32>> Lp,Rp;
	
	for(s32 i=l;i<mid;++i)
	{
		set_3d_VG_point(Lp,i,i+1,data,P1,Px,Py);
	}
	for(s32 i=r;i>mid;--i)
	{
		set_3d_VG_point(Rp,i,i-1,data,P1,Px,Py);
	}
	
	for(const auto &i:Lp) P1[mid].emplace_back(i);
	for(const auto &i:Rp) P1[mid].emplace_back(i);
	set_dis(P1[mid]);
	
	recursive_set_3d_VG_point(l,mid-1,data,P1,Px,Py);
	recursive_set_3d_VG_point(mid+1,r,data,P1,Px,Py);
}

void recursive_set_2D_VG(s32 pl,s32 pr,s32 sl, s32 sr,std::vector<std::pair<u32,u32>> &S,std::vector<Statemant_2D_VG> &state,std::vector<u32> &Px,s8 is_rev=0)
{
	if(pl>=pr)return;
	s32 pmid=(pl+pr)/2;
	
	s32 sl_mid,sr_mid;
	
	std::set<u32> ST;
	for(sl_mid=sl;state[sl_mid].a<Px[pmid];++sl_mid)
	{
		if(state[sl_mid].type==2)
		{
			ST.insert(state[sl_mid].b1);
		}
		else
		{
			auto it_l=ST.lower_bound(state[sl_mid].b1);
			auto it_r=ST.upper_bound(state[sl_mid].b2);
			while(it_l!=it_r)
			{
				auto tmp=it_l++;
				ST.erase(tmp);
			}
		}
	}
	if(is_rev)
	{
		for(auto i:ST)
		{
			S.emplace_back(i,Px[pmid]);
		}
	}
	else
	{
		for(auto i:ST)
		{
			S.emplace_back(Px[pmid],i);
		}
	}
	ST.clear();
	for(sr_mid=sr;state[sr_mid].a>Px[pmid];--sr_mid)
	{
		if(state[sr_mid].type==2)
		{
			ST.insert(state[sr_mid].b1);
		}
		else
		{
			auto it_l=ST.lower_bound(state[sr_mid].b1);
			auto it_r=ST.upper_bound(state[sr_mid].b2);
			while(it_l!=it_r)
			{
				auto tmp=it_l++;
				ST.erase(tmp);
			}
		}
	}
	if(is_rev)
	{
		for(auto i:ST)
		{
			S.emplace_back(i,Px[pmid]);
		}
	}
	else
	{
		for(auto i:ST)
		{
			S.emplace_back(Px[pmid],i);
		}
	}
	
	recursive_set_2D_VG(pl,pmid-1,sl,sl_mid-1,S,state,Px,is_rev);
	recursive_set_2D_VG(pmid+1,pr,sr_mid+1,sr,S,state,Px,is_rev);
}

inline void build_2D_VG_X(s32 lay,const DataSet &data,std::vector<std::pair<u32,u32>> *P1,std::vector<std::pair<u32,u32>> *P2,std::vector<s64> &Px,std::vector<s64> &Py)
{
	std::vector<Statemant_2D_VG> state;
	
	std::vector<u32> X;
	
	for(const auto &p:P1[lay])
	{
		state.emplace_back(2,p.first,p.second);
		X.emplace_back(p.first);
	}
	
	set_dis(X);
	
	for(const auto &o:data.Obstacles[lay])
	{
		s32 x1=get_dis(Px,o.first.x);
		s32 x2=s32(get_upper_dis(Px,o.second.x))-1;
		s32 y1=get_dis(Py,o.first.y);
		s32 y2=get_upper_dis(Py,o.second.y);
		
		if(x1<x2&&y1<y2)
		{
			state.emplace_back(3,x1,y1,y2);
			state.emplace_back(1,x2,y1,y2);
		}
	}
	
	sort(state.begin(),state.end());
	
	recursive_set_2D_VG(0,s32(X.size())-1,0,s32(state.size())-1,P2[lay],state,X);
}

inline void build_2D_VG_Y(s32 lay,const DataSet &data,std::vector<std::pair<u32,u32>> *P1,std::vector<std::pair<u32,u32>> *P2,std::vector<s64> &Px,std::vector<s64> &Py)
{
	std::vector<Statemant_2D_VG> state;
	
	std::vector<u32> Y;
	
	for(const auto &p:P1[lay])
	{
		state.emplace_back(2,p.second,p.first);
		Y.emplace_back(p.second);
	}
	
	set_dis(Y);
	
	for(const auto &o:data.Obstacles[lay])
	{
		s32 x1=get_dis(Px,o.first.x);
		s32 x2=get_upper_dis(Px,o.second.x);
		s32 y1=get_dis(Py,o.first.y);
		s32 y2=s32(get_upper_dis(Py,o.second.y))-1;
		
		if(x1<x2&&y1<y2)
		{
			state.emplace_back(3,y1,x1,x2);
			state.emplace_back(1,y2,x1,x2);
		}
	}
	
	sort(state.begin(),state.end());
	
	recursive_set_2D_VG(0,s32(Y.size())-1,0,s32(state.size())-1,P2[lay],state,Y,1);
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
	std::vector<std::pair<u32,u32>> P2[LIMIT_LAYER];
	
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		for(const auto &line:xLine[lay])
		{
			u32 x=get_dis(Px,std::get<0>(line));
			if(Px[x]<data.boundary.first.x||data.boundary.second.x<Px[x])
				continue;
			u32 l=get_dis(Py,std::get<1>(line));
			u32 r=get_dis(Py,std::get<2>(line));
			for(;l<=r;++l)
			{
				if(Py[l]<data.boundary.first.y||data.boundary.second.y<Py[l])
					continue;
				P1[lay].emplace_back(x,l);
			}
		}
		for(const auto &line:yLine[lay])
		{
			u32 y=get_dis(Py,std::get<0>(line));
			if(Py[y]<data.boundary.first.y||data.boundary.second.y<Py[y])
				continue;
			u32 l=get_dis(Px,std::get<1>(line));
			u32 r=get_dis(Px,std::get<2>(line));
			for(;l<=r;++l)
			{
				if(Px[l]<data.boundary.first.x||data.boundary.second.x<Px[l])
					continue;
				P1[lay].emplace_back(l,y);
			}
		}
		set_dis(P1[lay]);
		cout<<"P1["<<lay<<"].size(): "<<P1[lay].size()<<endl;
	}
	
	cout<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	recursive_set_3d_VG_point(1,data.metal_layers,data,P1,Px,Py);
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		cout<<"P1["<<lay<<"].size(): "<<P1[lay].size()<<endl;
	}
	cout<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		build_2D_VG_X(lay,data,P1,P2,Px,Py);
		cout<<"P2["<<lay<<"].size(): "<<P2[lay].size()<<endl;
		build_2D_VG_Y(lay,data,P1,P2,Px,Py);
		cout<<"  P2["<<lay<<"].size(): "<<P2[lay].size()<<endl;
		for(const auto &p:P2[lay])
		{
			P1[lay].emplace_back(p);
		}
		set_dis(P1[lay]);
	}
	
	
}