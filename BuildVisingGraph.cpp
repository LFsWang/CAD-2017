#include"DataLoader.h"
#include"BuildVisingGraph.h"
#include"SwapLineP1.h"
#include"BinaryIndexTree.h"
#include"DisjoinSet.h"
#include"allocator.h"
#include"showtime.h"
#include<omp.h>

const int debug_mode=0;

typedef std::set<u32,std::less<u32>,map_allocator<u32>> jinkela_set;

#include<future>
using std::ref;

template<typename T>
inline void set_dis(T &dis)
{
	std::sort(dis.begin(),dis.end());
	dis.resize(std::unique(dis.begin(),dis.end())-dis.begin());
}
template<typename T>
inline u32 get_dis(const std::vector<T> &dis,const T &val)
{
	return std::lower_bound(dis.begin(),dis.end(),val)-dis.begin();
}
template<typename T>
inline u32 get_upper_dis(const std::vector<T> &dis,const T &val)
{
	return std::upper_bound(dis.begin(),dis.end(),val)-dis.begin();
}

inline void get_original_PxPy(std::vector<s64> &Px,std::vector<s64> &Py,const DataSet &data)
{
	for(s32 i=1;i<=data.metal_layers;++i)
	{
		for(const auto &o:data.Obstacles[i])
		{
			Px.emplace_back(o.first.x);
			Px.emplace_back(o.second.x);
			Py.emplace_back(o.first.y);
			Py.emplace_back(o.second.y);
		}
		for(const auto &s:data.RoutedShape[i])
		{
			Px.emplace_back(s.first.x);
			Px.emplace_back(s.second.x);
			Py.emplace_back(s.first.y);
			Py.emplace_back(s.second.y);
		}
	}
	
	Px.emplace_back(data.boundary.first.x);
	Px.emplace_back(data.boundary.second.x);
	Py.emplace_back(data.boundary.first.y);
	Py.emplace_back(data.boundary.second.y);
	
	set_dis(Px);
	set_dis(Py);
}

inline void get_Line_swap_line(swape_line_P1 &SL,std::vector<Statemant_2D_VG> &line,std::vector<Statemant_2D_VG> &sline,std::vector<Statemant_2D_VG> &oline,std::vector<statementP1> &state,u32 l,u32 r,u32 al,u32 ar,u32 bl,u32 br)
{
	sort(state.begin(),state.end());
	
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
		s32 type=st.inout>0?3:1;
		for(const auto &seg:SL.segments)
		{
			line.emplace_back(type,st.a,seg.first,seg.second,st.seg_type);
		}
		if(SL.segments.size()&&al<=st.a&&st.a<=ar)
		{
			u32 L=SL.findL(SL.segments.front().first,l,r,1);
			u32 R=SL.findR(SL.segments.back().second,l,r,1);
			if(L<bl)L=bl;
			if(R>br)R=br;
			/*if(st.seg_type=='S')
			{
				if( L<SL.segments.front().first )
				{
					sline.emplace_back(1,st.a,L,SL.segments.front().first,'p');
					//oline.emplace_back(1,st.a,L,SL.segments.front().first,'p');
				}
				if( R>SL.segments.back().second )
				{
					sline.emplace_back(1,st.a,SL.segments.back().second,R,'p');
					//oline.emplace_back(1,st.a,SL.segments.back().second,R,'p');
				}
			}
			else
			{*/
				if( L<SL.segments.front().first )
				{
					oline.emplace_back(1,st.a,L,SL.segments.front().first,'p');
				}
				if( R>SL.segments.back().second )
				{
					oline.emplace_back(1,st.a,SL.segments.back().second,R,'p');
				}
			//}
		}
		SL.segments.clear();
	}
}

inline void get_xyLine_singal_layer(s32 lay,std::vector<Statemant_2D_VG> &xLine,std::vector<Statemant_2D_VG> &yLine,std::vector<Statemant_2D_VG> &sxLine,std::vector<Statemant_2D_VG> &syLine,std::vector<Statemant_2D_VG> &oxLine,std::vector<Statemant_2D_VG> &oyLine,const DataSet &data,const std::vector<s64> &Px,const std::vector<s64> &Py,swape_line_P1 &SLX,swape_line_P1 &SLY)
{
	std::vector<statementP1> Xstate,Ystate;
	
	for(const auto &o:data.Obstacles[lay])
	{
		u32 x1=get_dis(Px,o.first.x);
		u32 x2=get_dis(Px,o.second.x);
		u32 y1=get_dis(Py,o.first.y);
		u32 y2=get_dis(Py,o.second.y);
		
		Xstate.emplace_back(x1,y1,y2, 1,'O');
		Xstate.emplace_back(x2,y1,y2,-1,'O');
		Ystate.emplace_back(y1,x1,x2, 1,'O');
		Ystate.emplace_back(y2,x1,x2,-1,'O');
	}
	for(const auto &o:data.RoutedShape[lay])
	{
		u32 x1=get_dis(Px,o.first.x);
		u32 x2=get_dis(Px,o.second.x);
		u32 y1=get_dis(Py,o.first.y);
		u32 y2=get_dis(Py,o.second.y);
		
		Xstate.emplace_back(x1,y1,y2, 1,'S');
		Xstate.emplace_back(x2,y1,y2,-1,'S');
		Ystate.emplace_back(y1,x1,x2, 1,'S');
		Ystate.emplace_back(y2,x1,x2,-1,'S');
	}
	
	u32 lx = get_dis(Px,data.boundary.first.x);
	u32 ly = get_dis(Py,data.boundary.first.y);
	u32 ux = get_dis(Px,data.boundary.second.x);
	u32 uy = get_dis(Py,data.boundary.second.y);
	
	std::future<void> fx = std::async(get_Line_swap_line,ref(SLX),ref(xLine),ref(sxLine),ref(oxLine),ref(Xstate),0,Py.size()-1,lx,ux,ly,uy);
	std::future<void> fy = std::async(get_Line_swap_line,ref(SLY),ref(yLine),ref(syLine),ref(oyLine),ref(Ystate),0,Px.size()-1,ly,uy,lx,ux);
	
	fx.wait();
	fy.wait();
	
	//get_Line_swap_line(xLine,sxLine,oxLine,Xstate,0,Py.size()-1,lx,ux,ly,uy);
	//get_Line_swap_line(yLine,syLine,oyLine,Ystate,0,Px.size()-1,ly,uy,lx,ux);
}

inline void get_xyLine_singal_layer_future(std::vector<Statemant_2D_VG> *xLine,std::vector<Statemant_2D_VG> *yLine,std::vector<Statemant_2D_VG> *sxLine,std::vector<Statemant_2D_VG> *syLine,std::vector<Statemant_2D_VG> *oxLine,std::vector<Statemant_2D_VG> *oyLine,const DataSet &data,const std::vector<s64> &Px,const std::vector<s64> &Py)
{
	std::vector<std::future<void>> task;
	
	swape_line_P1 SLX[10+1];
	swape_line_P1 SLY[10+1];
	
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		SLX[lay].init(Py.size());
		SLY[lay].init(Px.size());
		task.emplace_back(std::async(get_xyLine_singal_layer,lay,ref(xLine[lay]),ref(yLine[lay]),ref(sxLine[lay]),ref(syLine[lay]),ref(oxLine[lay]),ref(oyLine[lay]),ref(data),ref(Px),ref(Py),ref(SLX[lay]),ref(SLY[lay])));
	}
	for(auto &f:task)
	{
		f.wait();
	}
}

inline void get_PxPy(std::vector<s64> &Px,std::vector<s64> &Py,const DataSet &data,std::vector<Statemant_2D_VG> *xLine,std::vector<Statemant_2D_VG> *yLine,std::vector<Statemant_2D_VG> *sxLine,std::vector<Statemant_2D_VG> *syLine,std::vector<Statemant_2D_VG> *oxLine,std::vector<Statemant_2D_VG> *oyLine)
{
	
	std::vector<s64> nPx,nPy;
	
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		for(const auto &i:xLine[lay])
		{
			nPx.emplace_back(Px[i.a]);
			nPy.emplace_back(Py[i.b1]);
			nPy.emplace_back(Py[i.b2]);
		}
		
		for(const auto &i:yLine[lay])
		{
			nPy.emplace_back(Py[i.a]);
			nPx.emplace_back(Px[i.b1]);
			nPx.emplace_back(Px[i.b2]);
		}
		for(const auto &i:data.RoutedVia[lay])
		{
			nPx.emplace_back(i.x);
			nPy.emplace_back(i.y);
		}
	}
	nPx.emplace_back(data.boundary.first.x);
	nPx.emplace_back(data.boundary.second.x);
	nPy.emplace_back(data.boundary.first.y);
	nPy.emplace_back(data.boundary.second.y);
	
	set_dis(nPx);
	set_dis(nPy);
	
	#pragma omp parallel for num_threads(data.metal_layers)
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		for(auto &i:xLine[lay])
		{
			i.a=get_dis(nPx,Px[i.a]);
			i.b1=get_dis(nPy,Py[i.b1]);
			i.b2=get_dis(nPy,Py[i.b2]);
		}
		
		for(auto &i:yLine[lay])
		{
			i.a=get_dis(nPy,Py[i.a]);
			i.b1=get_dis(nPx,Px[i.b1]);
			i.b2=get_dis(nPx,Px[i.b2]);
		}
		
		for(auto &i:sxLine[lay])
		{
			i.a=get_dis(nPx,Px[i.a]);
			i.b1=get_dis(nPy,Py[i.b1]);
			i.b2=get_dis(nPy,Py[i.b2]);
		}
		
		for(auto &i:syLine[lay])
		{
			i.a=get_dis(nPy,Py[i.a]);
			i.b1=get_dis(nPx,Px[i.b1]);
			i.b2=get_dis(nPx,Px[i.b2]);
		}
		
		for(auto &i:oxLine[lay])
		{
			i.a=get_dis(nPx,Px[i.a]);
			i.b1=get_dis(nPy,Py[i.b1]);
			i.b2=get_dis(nPy,Py[i.b2]);
		}
		
		for(auto &i:oyLine[lay])
		{
			i.a=get_dis(nPy,Py[i.a]);
			i.b1=get_dis(nPx,Px[i.b1]);
			i.b2=get_dis(nPx,Px[i.b2]);
		}
	}
	
	Px.swap(nPx);
	Py.swap(nPy);
}

inline void get_original_P1(std::vector<std::pair<u32,u32>> *P1,s32 metal_layers,std::vector<Statemant_2D_VG> *xLine,std::vector<Statemant_2D_VG> *yLine,const std::vector<s64> &Px,const std::vector<s64> &Py,const DataSet &data)
{
	for(s32 lay=1;lay<=metal_layers;++lay)
	{
		for(const auto &line:xLine[lay])
		{
			u32 x=line.a;
			u32 l=line.b1;
			u32 r=line.b2;
			
			//std::cerr<<x<<" ("<<l<<","<<r<<")\n";
			//std::cerr<<"----"<<Px[x]<<" ("<<Py[l]<<","<<Py[r]<<") "<<line.type<<' '<<char(line.seg_type)<<endl;
			
			P1[lay].emplace_back(x,l);
			P1[lay].emplace_back(x,r);
		}
		for(const auto &line:yLine[lay])
		{
			u32 y=line.a;
			u32 l=line.b1;
			u32 r=line.b2;
			
			P1[lay].emplace_back(l,y);
			P1[lay].emplace_back(r,y);
		}
		
		for(const auto &p:data.RoutedVia[lay])
		{
			P1[lay].emplace_back(get_dis(Px,p.x),get_dis(Py,p.y));
			P1[lay+1].emplace_back(get_dis(Px,p.x),get_dis(Py,p.y));
		}
		
		set_dis(P1[lay]);
	}
}

inline void set_3d_VG_point(std::vector<std::pair<u32,u32>> &S,s32 la1,s32 la2,const DataSet &data,std::vector<std::pair<u32,u32>> *P1,std::vector<s64> &Px,std::vector<s64> &Py,std::vector<std::pair<u32,u32>> *P2=nullptr)
{
	std::vector<Statemant_2D_VG> state;
	
	for(const auto &i:P1[la1]) S.emplace_back(i);
	
	state.reserve(S.size()+data.Obstacles[la2].size()*2);
	
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
		
		if(o.first.x!=o.second.x&&x1<=x2&&y1<y2)
		{
			state.emplace_back(3,x1,y1,y2);
			state.emplace_back(1,x2,y1,y2);
		}
	}
	
	sort(state.begin(),state.end());
	
	BIT ST;
	ST.init(Py.size()+2);
	
	std::vector<std::pair<u32,u32>> nS;
	nS.reserve(S.size());
	
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
	if(P2==nullptr) return;
	
	for(const auto &p:S)
	{
		P2[la2].emplace_back(p);
	}
}

//*
void project_point_on_all_layer(s32 l,s32 r,const DataSet &data,std::vector<std::pair<u32,u32>> *P1,std::vector<std::pair<u32,u32>> *P2,std::vector<s64> &Px,std::vector<s64> &Py)
{
	std::vector<std::pair<u32,u32>> Lp,Rp;
	
	for(s32 i=l;i<r;++i)
	{
		set_3d_VG_point(Lp,i,i+1,data,P1,Px,Py,P2);
	}
	for(s32 i=r;i>l;--i)
	{
		set_3d_VG_point(Rp,i,i-1,data,P1,Px,Py,P2);
	}
	for(s32 i=l;i<=r;++i)
	{
		for(const auto &j:P2[i])
		{
			P1[i].emplace_back(j);
		}
		P2[i].clear();
	}
}
//*/
	
void one_way_point_project_swap_line(const Statemant_2D_VG &st,jinkela_set &ST,std::vector<std::pair<u32,u32>> &S2,u32 L,u32 R,bool is_rev=0)
{
	if(st.type==2)
	{
		ST.emplace(st.b1);
	}
	else
	{
		
		if(st.a<L||R<st.a) return;
		
		auto it_l=ST.lower_bound(st.b1);
		auto it_r=ST.upper_bound(st.b2);
		
		if(!is_rev)
		{
			for(auto it=it_l;it!=it_r;++it)
			{
				S2.emplace_back(st.a,*it);
			}
		}
		else
		{
			for(auto it=it_l;it!=it_r;++it)
			{
				S2.emplace_back(*it,st.a);
			}
		}
		ST.erase(it_l,it_r);
	}
}

void single_layer_point_project(const std::vector<std::pair<u32,u32>> &S,std::vector<Statemant_2D_VG> &xLine,std::vector<Statemant_2D_VG> &yLine,std::vector<std::pair<u32,u32>> &S2,u32 x1,u32 x2,u32 y1,u32 y2)
{
	std::vector<Statemant_2D_VG> Xstate=xLine;
	std::vector<Statemant_2D_VG> Ystate=yLine;
	
	Xstate.reserve(xLine.size()+S.size());
	Ystate.reserve(yLine.size()+S.size());
	
	for(const auto &p:S)
	{
		Xstate.emplace_back(2,p.first,p.second);
		Ystate.emplace_back(2,p.second,p.first);
	}
	
	sort(Xstate.begin(),Xstate.end());
	sort(Ystate.begin(),Ystate.end());
	
	//std::cerr<<map_allocator_max<<endl;
	jinkela_set ST,RST;
	
	for(size_t i=0;i<Xstate.size();++i)
	{
		size_t ri=Xstate.size()-i-1;
		//if(Xstate[i].seg_type!='B'||Xstate[i].type!=3)
			one_way_point_project_swap_line(Xstate[i],ST,S2,x1,x2);
		//if(Xstate[ri].seg_type!='B'||Xstate[ri].type!=1)
			one_way_point_project_swap_line(Xstate[ri],RST,S2,x1,x2);
	}
	
	jinkela_set ST2,RST2;
	
	for(size_t i=0;i<Ystate.size();++i)
	{
		size_t ri=Ystate.size()-i-1;
		//if(Ystate[i].seg_type!='B'||Ystate[i].type!=3)
			one_way_point_project_swap_line(Ystate[i],ST2,S2,y1,y2,1);
		//if(Ystate[ri].seg_type!='B'||Ystate[ri].type!=1)
			one_way_point_project_swap_line(Ystate[ri],RST2,S2,y1,y2,1);
	}
	// to do thread
}

inline void point_project_to_XYLine_singal_layer(std::vector<std::pair<u32,u32>> &P1,std::vector<std::pair<u32,u32>> &P2,std::vector<std::pair<u32,u32>> &P3,std::vector<std::pair<u32,u32>> &P4,std::vector<Statemant_2D_VG> &xLine,std::vector<Statemant_2D_VG> &yLine,u32 x1,u32 x2,u32 y1,u32 y2)
{
	/*
	xLine.emplace_back(1,x1,y1,y2,'B');
	xLine.emplace_back(3,x2,y1,y2,'B');
	yLine.emplace_back(1,y1,x1,x2,'B');
	yLine.emplace_back(3,y2,x1,x2,'B');
	*/
	
	single_layer_point_project(P1,xLine,yLine,P2,x1,x2,y1,y2);
	
	for(const auto &p:P1)
	{
		P2.emplace_back(p);
	}
	
	for(const auto &p:P3)
	{
		P2.emplace_back(p);
	}
	std::vector<std::pair<u32,u32>>().swap(P3);
	for(const auto &p:P4)
	{
		P2.emplace_back(p);
	}
	std::vector<std::pair<u32,u32>>().swap(P4);
	
	set_dis(P2);
	P1.clear();
	P1.reserve(P2.size());
	for(const auto &p:P2)
	{
		if(x1<=p.first&&p.first<=x2&&y1<=p.second&&p.second<=y2)
			P1.emplace_back(p);
	}
	P2.clear();
}

inline void point_project_to_XYLine(std::vector<std::pair<u32,u32>> *P1,std::vector<std::pair<u32,u32>> *P2,std::vector<std::pair<u32,u32>> *P3,std::vector<std::pair<u32,u32>> *P4,const DataSet &data,std::vector<Statemant_2D_VG> *xLine,std::vector<Statemant_2D_VG> *yLine,const std::vector<s64> &Px,const std::vector<s64> &Py)
{
	u32 x1=get_dis(Px,data.boundary.first.x);
	u32 x2=get_dis(Px,data.boundary.second.x);
	u32 y1=get_dis(Py,data.boundary.first.y);
	u32 y2=get_dis(Py,data.boundary.second.y);
	
	map_allocator_max = 0;
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		map_allocator_max = std::max(P1[lay].size(),map_allocator_max);
	}
	map_allocator_max += 4;
	
	std::vector<std::future<void>> task;
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		P2[lay].reserve(P1[lay].size()*4+P3[lay].size()+P4[lay].size());
		task.emplace_back(std::async(point_project_to_XYLine_singal_layer,ref(P1[lay]),ref(P2[lay]),ref(P3[lay]),ref(P4[lay]),ref(xLine[lay]),ref(yLine[lay]),x1,x2,y1,y2));
	}
	for(auto &f:task)
	{
		f.wait();
	}
}

void SO_point_project_swap_line(const Statemant_2D_VG &st,jinkela_set &ST,std::vector<std::pair<u32,u32>> &S2,bool is_rev=0)
{
	if(st.type==2)
	{
		ST.emplace(st.b1);
	}
	else
	{
		auto it_l=ST.lower_bound(st.b1);
		auto it_r=ST.upper_bound(st.b2);
		
		if(st.seg_type=='p')
		{
			if(!is_rev)
			{
				for(auto it=it_l;it!=it_r;++it)
				{
					S2.emplace_back(st.a,*it);
				}
			}
			else
			{
				for(auto it=it_l;it!=it_r;++it)
				{
					S2.emplace_back(*it,st.a);
				}
			}
		}
		ST.erase(it_l,it_r);
	}
}

inline void point_project_to_soxyLine_singal_layer(const std::vector<std::pair<u32,u32>> &S,std::vector<Statemant_2D_VG> &xLine,std::vector<Statemant_2D_VG> &yLine,std::vector<Statemant_2D_VG> &pxLine,std::vector<Statemant_2D_VG> &pyLine,std::vector<std::pair<u32,u32>> &S2)
{
	std::vector<Statemant_2D_VG> LXstate=xLine;
	std::vector<Statemant_2D_VG> RXstate=xLine;
	std::vector<Statemant_2D_VG> LYstate=yLine;
	std::vector<Statemant_2D_VG> RYstate=yLine;
	
	LXstate.reserve(xLine.size()+pxLine.size()+S.size());
	RXstate.reserve(xLine.size()+pxLine.size()+S.size());
	LYstate.reserve(yLine.size()+pyLine.size()+S.size());
	RYstate.reserve(yLine.size()+pyLine.size()+S.size());
	
	for(const auto &st:pxLine)
	{
		RXstate.emplace_back(st);
		LXstate.emplace_back(st);
		LXstate.back().type=3;
	}
	
	for(const auto &st:pyLine)
	{
		RYstate.emplace_back(st);
		LYstate.emplace_back(st);
		LYstate.back().type=3;
	}
	
	for(const auto &p:S)
	{
		LXstate.emplace_back(2,p.first,p.second);
		RXstate.emplace_back(2,p.first,p.second);
		LYstate.emplace_back(2,p.second,p.first);
		RYstate.emplace_back(2,p.second,p.first);
	}
	
	sort(LXstate.begin(),LXstate.end());
	sort(RXstate.begin(),RXstate.end());
	sort(LYstate.begin(),LYstate.end());
	sort(RYstate.begin(),RYstate.end());
	
	jinkela_set XST,XRST,YST,YRST;
	for(const auto &st:LXstate)
	{
		SO_point_project_swap_line(st,XST,S2);
		//one_way_point_project_swap_line(st,XST,S2,0,2147483647);
	}
	
	for(auto it=RXstate.rbegin();it!=RXstate.rend();++it)
	{
		SO_point_project_swap_line(*it,XRST,S2);
		//one_way_point_project_swap_line(*it,XRST,S2,0,2147483647);
	}
	
	for(const auto &st:LYstate)
	{
		SO_point_project_swap_line(st,YST,S2,1);
		//one_way_point_project_swap_line(st,YST,S2,0,2147483647,1);
	}
	
	for(auto it=RYstate.rbegin();it!=RYstate.rend();++it)
	{
		SO_point_project_swap_line(*it,YRST,S2,1);
		//one_way_point_project_swap_line(*it,YRST,S2,0,2147483647,1);
	}
}

inline void point_project_to_soxyLine(std::vector<std::pair<u32,u32>> *P1,std::vector<std::pair<u32,u32>> *P2,const DataSet &data,std::vector<Statemant_2D_VG> *xLine,std::vector<Statemant_2D_VG> *yLine,std::vector<Statemant_2D_VG> *pxLine,std::vector<Statemant_2D_VG> *pyLine)
{
	map_allocator_max = 0;
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		map_allocator_max = std::max(map_allocator_max,P1[lay].size());
	}
	map_allocator_max += 4;
	
	std::vector<std::future<void>> task;
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		P2[lay].reserve(P1[lay].size()*4);
		task.emplace_back(std::async(point_project_to_soxyLine_singal_layer,ref(P1[lay]),ref(xLine[lay]),ref(yLine[lay]),ref(pxLine[lay]),ref(pyLine[lay]),ref(P2[lay])));
	}
	for(auto &f:task)
	{
		f.wait();
	}
}

inline void set_Pv(const std::vector<point3D> &V_set,s32 lay,std::vector<std::pair<u32,u32>> *P1,std::vector<size_t> *Pv)
{
	for(const auto &p:P1[lay])
	{
		Pv[lay].emplace_back(get_dis(V_set,point3D(p.first,p.second,lay)));
	}
	std::vector<std::pair<u32,u32>>().swap(P1[lay]);
}

inline void put_the_point_number(s32 metal_layers,std::vector<point3D> &V_set,std::vector<std::pair<u32,u32>> *P1,std::vector<size_t> *Pv)
{
	size_t cnt = 0;
	for(s32 lay=1;lay<=metal_layers;++lay)
	{
		cnt+=P1[lay].size();
	}
	
	V_set.clear();
	V_set.reserve(cnt);
	
	for(s32 lay=1;lay<=metal_layers;++lay)
	{
		//P1[lay].shrink_to_fit();
		for(const auto &p:P1[lay])
		{
			V_set.emplace_back(p.first,p.second,lay);
		}
	}
	set_dis(V_set);
	std::vector<std::future<void>> task;
	for(s32 lay=1;lay<=metal_layers;++lay)
	{
		Pv[lay].reserve(P1[lay].size());
		task.emplace_back(std::async(set_Pv,ref(V_set),lay,P1,Pv));
	}
	for(auto &f:task)
	{
		f.wait();
	}
}

inline void merge_same_via_point(DisjoinSet &DST,const DataSet &data,const std::vector<point3D> &V_set,const std::vector<s64> &Px,const std::vector<s64> &Py)
{
	for(s32 lay=1;lay<data.metal_layers;++lay)
	{
		for(const auto &p:data.RoutedVia[lay])
		{
			u32 x=get_dis(Px,p.x);
			u32 y=get_dis(Py,p.y);
			
			point3D P3D1(x,y,lay);
			point3D P3D2(x,y,lay+1);
			
			size_t p1=get_dis(V_set,P3D1);
			size_t p2=get_dis(V_set,P3D2);
			
			if(debug_mode)
			{
				if(!(p1<V_set.size()&&V_set[p1]==P3D1))std::cerr<<"Error V_set[p1]!=P3D1\n";
				if(!(p2<V_set.size()&&V_set[p2]==P3D2))std::cerr<<"Error V_set[p2]!=P3D2\n";
			}
			
			DST.U(p1,p2);
		}
	}
}

inline void set_shape_point_singal_layer(std::vector<bool> &ori_is_P,s32 lay,std::vector<Statemant_2D_VG> *xLine,std::vector<Statemant_2D_VG> *yLine,const std::vector<point3D> &V_set,DisjoinSet &DST)
{
	for(const auto &o:xLine[lay])
	{
		if(o.seg_type!='S')continue;
		u32 x1=o.a;
		u32 x2=o.a;
		u32 y1=o.b1;
		u32 y2=o.b2;
		
		point3D P3D1(x1,y1,lay);
		point3D P3D2(x2,y2,lay);
		
		size_t p1=get_dis(V_set,P3D1);
		size_t p2=get_dis(V_set,P3D2);
		
		DST.U(p1,p2);
		
		if(debug_mode)
		{
			if(!(p1<V_set.size()&&V_set[p1]==P3D1))std::cerr<<"Error V_set[p1]!=P3D1\n";
			if(!(p2<V_set.size()&&V_set[p2]==P3D2))std::cerr<<"Error V_set[p2]!=P3D2\n";
		}
		
		ori_is_P[p1]=true;
		ori_is_P[p2]=true;
	}
	for(const auto &o:yLine[lay])
	{
		if(o.seg_type!='S')continue;
		u32 y1=o.a;
		u32 y2=o.a;
		u32 x1=o.b1;
		u32 x2=o.b2;
		
		point3D P3D1(x1,y1,lay);
		point3D P3D2(x2,y2,lay);
		
		size_t p1=get_dis(V_set,P3D1);
		size_t p2=get_dis(V_set,P3D2);
		
		DST.U(p1,p2);
		
		
		if(debug_mode)
		{
			if(!(p1<V_set.size()&&V_set[p1]==P3D1))std::cerr<<"Error V_set[p1]!=P3D1\n";
			if(!(p2<V_set.size()&&V_set[p2]==P3D2))std::cerr<<"Error V_set[p2]!=P3D2\n";
		}
		
		ori_is_P[p1]=true;
		ori_is_P[p2]=true;
	}
}

inline void set_shape_point(std::vector<bool> &ori_is_P,s32 metal_layers,std::vector<Statemant_2D_VG> *xLine,std::vector<Statemant_2D_VG> *yLine,const std::vector<point3D> &V_set,DisjoinSet &DST)
{
	std::vector<std::future<void>> task;
	for(s32 lay=1;lay<=metal_layers;++lay)
	{
		task.emplace_back(std::async(set_shape_point_singal_layer,ref(ori_is_P),lay,xLine,yLine,ref(V_set),ref(DST)));
		//set_shape_point_singal_layer(ori_is_P,lay,xLine,yLine,V_set,DST);
	}
	for(auto &t:task)
	{
		t.wait();
	}
}

inline void merge_shape_point_swap_line(DisjoinSet &DST,std::vector<Statemant_2D_VG> &state,u32 bit_size)
{
	sort(state.begin(),state.end());
	
	std::map<u32,u32,std::less<u32>,map_allocator<std::pair<u32,u32>>> ST;
	
	BIT BIST;
	BIST.init(bit_size+2);
	
	for(const auto &st:state)
	{
		if(st.type==2)
		{
			if(BIST.get_sum(st.b1+1)||BIST.get_sum(st.b1))
			{
				auto it=ST.find(st.b1);
				if(it!=ST.end())
				{
					DST.U(it->second,st.b2);
					it->second = st.b2;
				}
				else
				{
					ST.emplace(st.b1,st.b2);
				}
			}
		}
		else
		{
			if(st.type==1)
			{
				BIST.add(st.b1+1,1);
				BIST.add(st.b2+1,-1);
			}
			else
			{
				BIST.add(st.b1+1,-1);
				BIST.add(st.b2+1,1);
			}
			if(st.type==1)continue;
			
			auto it_l=ST.lower_bound(st.b1);
			auto it_r=ST.upper_bound(st.b2);
			size_t p1 = it_l->second;
			
			for(auto it=it_l;it!=it_r;++it)
			{
				DST.U(p1,it->second);
			}
			ST.erase(it_l,it_r);
		}
	}
}

inline void merge_same_shape_point_singal_layer(DisjoinSet &DST,std::vector<Statemant_2D_VG> *xLine,std::vector<Statemant_2D_VG> *yLine,std::vector<size_t> *Pv,const std::vector<point3D> &V_set,const std::vector<s64> &Px,const std::vector<s64> &Py,s32 lay)
{
	std::vector<Statemant_2D_VG> Xstate;//global
	std::vector<Statemant_2D_VG> Ystate;
	
	Xstate.reserve(xLine[lay].size()+Pv[lay].size());
	Ystate.reserve(yLine[lay].size()+Pv[lay].size());
	
	for(const auto &st:xLine[lay])
	{
		if(st.seg_type!='S') continue;
		s32 type=st.type==3?1:3;
		Xstate.emplace_back(type,st.a,st.b1,st.b2,'S');
	}
	for(const auto &st:yLine[lay])
	{
		if(st.seg_type!='S') continue;
		s32 type=st.type==3?1:3;
		Ystate.emplace_back(type,st.a,st.b1,st.b2,'S');
	}
	
	for(auto p:Pv[lay])
	{
		Xstate.emplace_back(2,V_set[p].x,V_set[p].y,p);
		Ystate.emplace_back(2,V_set[p].y,V_set[p].x,p);
	}
	
	merge_shape_point_swap_line(DST,Xstate,Py.size());
	merge_shape_point_swap_line(DST,Ystate,Px.size());
}

inline void merge_same_shape_point(DisjoinSet &DST,std::vector<Statemant_2D_VG> *xLine,std::vector<Statemant_2D_VG> *yLine,std::vector<size_t> *Pv,const std::vector<point3D> &V_set,const std::vector<s64> &Px,const std::vector<s64> &Py,s32 metal_layers)
{
	map_allocator_max = 0;
	for(s32 lay=1;lay<=metal_layers;++lay)
	{
		map_allocator_max = std::max(Pv[lay].size(),map_allocator_max);
	}
	map_allocator_max += 4;
	
	std::vector<std::future<void>> task;
	for(s32 lay=1;lay<=metal_layers;++lay)
	{
		task.emplace_back(std::async(merge_same_shape_point_singal_layer,ref(DST),xLine,yLine,Pv,ref(V_set),ref(Px),ref(Py),lay));
		//merge_same_shape_point_singal_layer(DST,xLine,yLine,V,V_set,Px,Py,lay);
	}
	for(auto &t:task)
	{
		t.wait();
	}
}

inline size_t shrink_point(std::vector<size_t> &shrink_from,std::vector<bool> &is_pinv,std::vector<bool> &ori_is_P,DisjoinSet &DST,size_t sz)
{
	shrink_from.resize(sz);
	
	for(size_t i=0;i<sz;++i)
	{
		shrink_from[i]=DST.find(i);
	}
	
	std::vector<size_t> tmp=shrink_from;
	set_dis(tmp);
	
	is_pinv.clear();
	is_pinv.resize(tmp.size());
	
	#pragma omp parallel for num_threads(8)
	for(size_t i=0;i<sz;++i)
	{
		shrink_from[i]=get_dis(tmp,shrink_from[i]);
	}
		
	for(size_t i=0;i<sz;++i)
	{
		is_pinv[shrink_from[i]] = is_pinv[shrink_from[i]] || ori_is_P[i];
	}
	
	return tmp.size();
}

inline void build_2D_edge_swape_line(const std::vector<size_t> &shrink_from,std::vector<Statemant_2D_VG> &state,std::vector<Edge> &edge,bool is_rev=0)
{
	sort(state.begin(),state.end());
	
	std::map<u32,u32,std::less<u32>,map_allocator<std::pair<u32,u32>>> ST;
	
	u8 edge_type = is_rev ? 'V' : 'H';
	
	for(const auto &st:state)
	{
		if(st.type==2)
		{
			auto it=ST.find(st.b1);
			if(it!=ST.end())
			{
				if(shrink_from[st.b2]!=shrink_from[it->second])
				{
					edge.emplace_back(st.b2,it->second,edge_type);
					edge.emplace_back(it->second,st.b2,edge_type);
				}
				it->second = st.b2;
			}
			else
			{
				ST.emplace(st.b1,st.b2);
			}
		}
		else
		{
			auto it_l=ST.lower_bound(st.b1);
			auto it_r=ST.upper_bound(st.b2);
			ST.erase(it_l,it_r);
		}
	}
	
}

inline void build_2D_edge_singal_layer(s32 lay,const DataSet &data,const std::vector<size_t> &shrink_from,std::vector<size_t> *Pv,std::vector<s64> &Px,std::vector<s64> &Py,std::vector<Edge> &edgeX,std::vector<Edge> &edgeY,const std::vector<point3D> &V_set,std::vector<Statemant_2D_VG> &Xstate,std::vector<Statemant_2D_VG> &Ystate)
{
	
	for(auto p:Pv[lay])
	{
		Xstate.emplace_back(2,V_set[p].x,V_set[p].y,p);
		Ystate.emplace_back(2,V_set[p].y,V_set[p].x,p);
	}

	for(const auto &o:data.Obstacles[lay])
	{
		s32 x1=get_upper_dis(Px,o.first.x);
		s32 x2=s32(get_upper_dis(Px,o.second.x))-1;
		s32 y1=get_upper_dis(Py,o.first.y);
		s32 y2=s32(get_upper_dis(Py,o.second.y))-1;
		
		if(x2>=0&&Px[x2]==o.second.x)--x2;
		if(y2>=0&&Py[y2]==o.second.y)--y2;
		
		if(x1<=x2&&y1<=y2)
		{
			Xstate.emplace_back(3,x1,y1,y2);
			Xstate.emplace_back(1,x2,y1,y2);
			
			Ystate.emplace_back(3,y1,x1,x2);
			Ystate.emplace_back(1,y2,x1,x2);
		}
	}
	
	std::future<void> XF(std::async(build_2D_edge_swape_line,ref(shrink_from),ref(Xstate),ref(edgeX),0));
	std::future<void> YF(std::async(build_2D_edge_swape_line,ref(shrink_from),ref(Ystate),ref(edgeY),1));
	
	XF.wait();
	YF.wait();
	
	//build_2D_edge_swape_line(shrink_from,Xstate,edgeX);
	//build_2D_edge_swape_line(shrink_from,Ystate,edgeY,1);
}

inline void build_2D_edge(const DataSet &data,const std::vector<size_t> &shrink_from,std::vector<size_t> *Pv,std::vector<s64> &Px,std::vector<s64> &Py,std::vector<Edge> &edge,const std::vector<point3D> &V_set)
{
	map_allocator_max = 0;
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		map_allocator_max = std::max(Pv[lay].size(),map_allocator_max);
	}
	map_allocator_max += 4;
	
	std::vector<std::future<void>> task;
	
	std::vector<Edge> edgeX[10+1],edgeY[10+1];
	std::vector<Statemant_2D_VG> Xstate[10+1];
	std::vector<Statemant_2D_VG> Ystate[10+1];
	
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		Xstate[lay].reserve(data.Obstacles[lay].size()*2+Pv[lay].size());
		Ystate[lay].reserve(data.Obstacles[lay].size()*2+Pv[lay].size());
		
		edgeX[lay].reserve(Pv[lay].size()*4);
		edgeY[lay].reserve(Pv[lay].size()*4);
		
		task.emplace_back(std::async(build_2D_edge_singal_layer,lay,ref(data),ref(shrink_from),Pv,ref(Px),ref(Py),ref(edgeX[lay]),ref(edgeY[lay]),ref(V_set),ref(Xstate[lay]),ref(Ystate[lay])));
		//build_2D_edge_singal_layer(lay,data,shrink_from,V,Px,Py,edgeX[lay],edgeY[lay],V_set);
	}
	
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		task[lay-1].wait();
		for(auto e:edgeX[lay])
		{
			edge.emplace_back(e);
		}
		edgeX[lay]=std::vector<Edge>();
		for(auto e:edgeY[lay])
		{
			edge.emplace_back(e);
		}
		edgeY[lay]=std::vector<Edge>();
	}
}

inline void build_via_edge_swap_line(s32 la1,std::vector<Statemant_2D_VG> &state_la2,std::vector<size_t> *Pv,std::vector<size_t> &res,u32 bit_size,const std::vector<point3D> &V_set)
{
	for(const auto p:res)
	{
		state_la2.emplace_back(2,V_set[p].x,V_set[p].y,p);
	}
	
	res.clear();
	
	for(const auto p:Pv[la1])
	{
		state_la2.emplace_back(2,V_set[p].x,V_set[p].y,p);
	}
	
	sort(state_la2.begin(),state_la2.end());
	
	BIT ST;
	ST.init(bit_size+2);
	
	for(const auto &st:state_la2)
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
					res.emplace_back(st.b2);
				}
				break;
			}
			case 3:{
				ST.add(st.b1+1,1);
				ST.add(st.b2+1,-1);
				break;
			}
			default:
				std::cerr<<"Error!\n";
		}
	}
}

inline void build_via_edge(const DataSet &data,const std::vector<size_t> &shrink_from,std::vector<size_t> *Pv,std::vector<s64> &Px,std::vector<s64> &Py,std::vector<Edge> &edge,const std::vector<point3D> &V_set)
{
	std::vector<size_t> res[2];
	std::vector<Statemant_2D_VG> state;
	
	for(s32 lay=2;lay<=data.metal_layers;++lay)
	{
		state.clear();
		for(const auto &o:data.Obstacles[lay])
		{
			s32 x1=get_dis(Px,o.first.x);
			s32 x2=s32(get_upper_dis(Px,o.second.x))-1;
			s32 y1=get_dis(Py,o.first.y);
			s32 y2=get_upper_dis(Py,o.second.y);
			
			if(x1<=x2&&y1<y2)
			{
				state.emplace_back(3,x1,y1,y2);
				state.emplace_back(1,x2,y1,y2);
			}
		}
		
		build_via_edge_swap_line(lay-1,state,Pv,res[lay&1],Py.size(),V_set);
		res[(lay+1)&1].clear();
		
		for(const auto p:res[lay&1])
		{
			auto P3D2=point3D(V_set[p].x,V_set[p].y,lay);
			size_t p2=get_dis(V_set,P3D2);
			if(V_set[p2]==P3D2)
			{
				if(shrink_from[p]!=shrink_from[p2])
				{
					edge.emplace_back(p,p2,'Z');
					edge.emplace_back(p2,p,'Z');
				}
			}
			else
			{
				res[(lay+1)&1].emplace_back(p);
			}
		}
	}
}

inline void set_edge(size_t l,size_t r,s64 viacost,std::vector<Edge> &edge,std::vector<s64> &Px,std::vector<s64> &Py,const std::vector<point3D> &V_set)
{
	for(;l<r;++l)
	{
		auto &e=edge[l];
		if(e.type=='Z')
		{
			s64 z1=V_set[e.ori_u].layer;
			s64 z2=V_set[e.ori_v].layer;
			e.cost=viacost*std::abs(z1-z2);
		}
		else
		{
			s64 x1=Px[V_set[e.ori_u].x];
			s64 x2=Px[V_set[e.ori_v].x];
			s64 y1=Py[V_set[e.ori_u].y];
			s64 y2=Py[V_set[e.ori_v].y];
			e.cost=std::abs(x1-x2)+std::abs(y1-y2);
		}
	}
}

inline void set_edge_and_graph(s64 viacost,size_t N,std::vector<std::vector<size_t>> &G,std::vector<Edge> &edge,const std::vector<size_t> &shrink_from,std::vector<s64> &Px,std::vector<s64> &Py,const std::vector<point3D> &V_set)
{
	size_t edge_div = edge.size()/7, edge_cnt = 0;
	std::vector<std::future<void>> task;
	for(s32 i=0;i<7;++i)
	{
		task.emplace_back(std::async(set_edge,edge_cnt,i==6?(edge.size()):(edge_cnt+edge_div),viacost,ref(edge),ref(Px),ref(Py),ref(V_set)));
		edge_cnt+=edge_div;
		if(debug_mode) std::cerr<<"edge_cnt: "<<edge_cnt<<'\n';
	}
	
	G.clear();
	G.resize(N);
	for(auto &g:G)
	{
		g.reserve(6);
	}
	
	for(size_t i=0;i<edge.size();++i)
	{
		auto &e=edge[i];
		
		e.u=shrink_from[e.ori_u];
		e.v=shrink_from[e.ori_v];
		
		if(debug_mode)
		{
			if(e.u==e.v)std::cerr<<"Error!\n";
		}
		
		G[e.u].emplace_back(i);
		
	}
	
	for(auto &t:task)
	{
		t.wait();
	}
}


void VisingGraph::build(const DataSet &data,bool is_not_connect=0)
{	
	static const int LIMIT_LAYER = DataSet::LIMIT_LAYER;
	
	std::vector<Statemant_2D_VG> sxLine[LIMIT_LAYER];
	std::vector<Statemant_2D_VG> oxLine[LIMIT_LAYER];
	std::vector<Statemant_2D_VG> syLine[LIMIT_LAYER];
	std::vector<Statemant_2D_VG> oyLine[LIMIT_LAYER];
	
	showclock("start");
	
	get_original_PxPy(Px,Py,data);
	
	showclock("get_original_PxPy");
	
	get_xyLine_singal_layer_future(xLine,yLine,sxLine,syLine,oxLine,oyLine,data,Px,Py);
	
	showclock("get_xyLine_singal_layer_future");
	
	get_PxPy(Px,Py,data,xLine,yLine,sxLine,syLine,oxLine,oyLine);
	showclock("get_PxPy");
	
	std::vector<std::pair<u32,u32>> P1[LIMIT_LAYER];
	std::vector<std::pair<u32,u32>> P2[LIMIT_LAYER];
	std::vector<std::pair<u32,u32>> P3[LIMIT_LAYER];
	std::vector<std::pair<u32,u32>> P4[LIMIT_LAYER];
	
	get_original_P1(P1,data.metal_layers,xLine,yLine,Px,Py,data);
	showclock("get_original_P1");
	
	//*
	if(debug_mode)
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::cerr<<"P1["<<lay<<"].size(): "<<P1[lay].size()<<endl;
	}
	if(debug_mode) std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	//*/
	
	if(is_not_connect){
		//build_2D_VG_point(data,P1,P2,Px,Py);
		//showclock("build_2D_VG_point");
	}
	
	point_project_to_XYLine(P1,P2,P3,P4,data,xLine,yLine,Px,Py);
	showclock("point_project_to_XYLine");
	
	project_point_on_all_layer(1,data.metal_layers,data,P1,P2,Px,Py);// insert more point
	showclock("project_point_on_all_layer");
	
	//*
	if(debug_mode)
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::cerr<<"P1["<<lay<<"].size(): "<<P1[lay].size()<<endl;
	}
	if(debug_mode) std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	//*/
	
	//*
	/*point_project_to_soxyLine(P1,P3,data,xLine,yLine,sxLine,syLine);
	showclock("shape: point_project_to_soxyLine");
	
	if(debug_mode)
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::cerr<<"P3["<<lay<<"].size(): "<<P3[lay].size()<<endl;
	}
	if(debug_mode) std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";*/
	
	point_project_to_soxyLine(P1,P4,data,xLine,yLine,oxLine,oyLine);
	showclock("obstacle: point_project_to_soxyLine");
	
	if(debug_mode)
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::cerr<<"P4["<<lay<<"].size(): "<<P4[lay].size()<<endl;
	}
	if(debug_mode) std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	//*/
	
	point_project_to_XYLine(P1,P2,P3,P4,data,xLine,yLine,Px,Py);
	if(debug_mode) showclock("point_project_to_XYLine");
	
	
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::vector<Statemant_2D_VG>().swap(sxLine[lay]);
		std::vector<Statemant_2D_VG>().swap(oxLine[lay]);
		std::vector<Statemant_2D_VG>().swap(syLine[lay]);
		std::vector<Statemant_2D_VG>().swap(oyLine[lay]);
	}
	
	//*
	if(debug_mode)
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::cerr<<"P1["<<lay<<"].size(): "<<P1[lay].size()<<endl;
	}
	if(debug_mode) std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	//*/
	
	std::vector<size_t> Pv[LIMIT_LAYER];
	
	put_the_point_number(data.metal_layers,V_set,P1,Pv);
	showclock("put_the_point_number");
	
	DST.init(V_set.size());
	
	std::vector<bool> ori_is_P(V_set.size());
	
	merge_same_via_point(DST,data,V_set,Px,Py);
	showclock("merge_same_via_point");
	
	set_shape_point(ori_is_P,data.metal_layers,xLine,yLine,V_set,DST);
	showclock("set_shape_point");
	
	merge_same_shape_point(DST,xLine,yLine,Pv,V_set,Px,Py,data.metal_layers);
	showclock("merge_same_shape_point");
	
	N = shrink_point(shrink_from,is_pinv,ori_is_P,DST,V_set.size());
	showclock("shrink_point");
	
	std::cerr<<N<<' '<<V_set.size()<<endl;
	
	edge.clear();
	edge.reserve(V_set.size()*12);
	showclock("edge.reserve(V_set.size()*12)");
	
	build_via_edge(data,shrink_from,Pv,Px,Py,edge,V_set);
	showclock("build_via_edge");
	
	build_2D_edge(data,shrink_from,Pv,Px,Py,edge,V_set);
	showclock("build_2D_edge");
	
	if(debug_mode)
	std::cerr<<"edge size: "<<edge.size()<<'\n';
	
	set_edge_and_graph(data.viacost,N,G,edge,shrink_from,Px,Py,V_set);
	showclock("set_edge_and_graph");
}

void VisingGraph::print_select_edges(const std::vector<std::size_t> &res,std::ofstream &fout)
{
	for(auto it:res)
	{
        auto i=it%2?it^1:it;
		if(edge[i].type=='Z'){
			for(auto lay=V_set[edge[i].ori_u].layer;lay<V_set[edge[i].ori_v].layer;++lay)
			fout<<"Via V"<<lay<<" ("<<Px[V_set[edge[i].ori_u].x]<<","<<Py[V_set[edge[i].ori_u].y]<<")\n";
		}
		else{
			fout<<edge[i].type<<"-line M"<<V_set[edge[i].ori_u].layer<<" ("<<Px[V_set[edge[i].ori_u].x]<<","<<Py[V_set[edge[i].ori_u].y]<<") ("<<Px[V_set[edge[i].ori_v].x]<<","<<Py[V_set[edge[i].ori_v].y]<<")\n";
		}
	}
}

// beta function

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
	
	std::future<void> L(std::async(recursive_set_3d_VG_point,l,mid-1,ref(data),P1,ref(Px),ref(Py)));
	std::future<void> R(std::async(recursive_set_3d_VG_point,mid+1,r,ref(data),P1,ref(Px),ref(Py)));
	
	L.wait();
	R.wait();
}

//*
void project_point_on_all_layer_beta(s32 l,s32 r,const DataSet &data,std::vector<std::pair<u32,u32>> *P1,std::vector<std::pair<u32,u32>> *P2,std::vector<s64> &Px,std::vector<s64> &Py)
{
	std::vector<std::pair<u32,u32>> Lp,Rp;
	
	for(s32 i=l;i<r;++i)
	{
		set_3d_VG_point(Lp,i,i+1,data,P1,Px,Py,P2);
	}
	for(s32 i=r;i>l;--i)
	{
		set_3d_VG_point(Rp,i,i-1,data,P1,Px,Py,P2);
	}
}
//*/

void recursive_set_2D_VG(s32 pl,s32 pr,s32 sl, s32 sr,std::vector<std::pair<u32,u32>> &S,std::vector<Statemant_2D_VG> &state,std::vector<u32> &Px,s8 is_rev)
{
	if(pl>=pr)return;
	s32 pmid=(pl+pr)/2;
	
	s32 sl_mid,sr_mid;
	
	std::set<u32> ST;
	for(sl_mid=sl;state[sl_mid].a<Px[pmid];++sl_mid)
	{
		if(state[sl_mid].type==2)
		{
			ST.emplace(state[sl_mid].b1);
		}
		else
		{
			auto it_l=ST.lower_bound(state[sl_mid].b1);
			auto it_r=ST.upper_bound(state[sl_mid].b2);
			ST.erase(it_l,it_r);
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
			ST.emplace(state[sr_mid].b1);
		}
		else
		{
			auto it_l=ST.lower_bound(state[sr_mid].b1);
			auto it_r=ST.upper_bound(state[sr_mid].b2);
			ST.erase(it_l,it_r);
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
inline void build_2D_VG_X_point(s32 lay,const DataSet &data,std::vector<std::pair<u32,u32>> *P1,std::vector<std::pair<u32,u32>> *P2,std::vector<s64> &Px,std::vector<s64> &Py)
{
	std::vector<Statemant_2D_VG> state;
	state.reserve(P1[lay].size()+data.Obstacles[lay].size()*2);
	
	std::vector<u32> X;
	X.reserve(P1[lay].size());
	
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
		s32 y2=s32(get_upper_dis(Py,o.second.y))-1;
		
		if(x1<=x2&&y1<=y2)
		{
			state.emplace_back(3,x1,y1,y2);
			state.emplace_back(1,x2,y1,y2);
		}
	}
	
	sort(state.begin(),state.end());
	
	recursive_set_2D_VG(0,s32(X.size())-1,0,s32(state.size())-1,P2[lay],state,X,0);
}
inline void build_2D_VG_Y_point(s32 lay,const DataSet &data,std::vector<std::pair<u32,u32>> *P1,std::vector<std::pair<u32,u32>> *P2,std::vector<s64> &Px,std::vector<s64> &Py)
{
	std::vector<Statemant_2D_VG> state;
	state.reserve(P1[lay].size()+data.Obstacles[lay].size()*2);
	
	std::vector<u32> Y;
	Y.reserve(P1[lay].size());
	
	for(const auto &p:P1[lay])
	{
		state.emplace_back(2,p.second,p.first);
		Y.emplace_back(p.second);
	}
	
	set_dis(Y);
	
	for(const auto &o:data.Obstacles[lay])
	{
		s32 x1=get_dis(Px,o.first.x);
		s32 x2=s32(get_upper_dis(Px,o.second.x))-1;
		s32 y1=get_dis(Py,o.first.y);
		s32 y2=s32(get_upper_dis(Py,o.second.y))-1;
		
		if(x1<=x2&&y1<=y2)
		{
			state.emplace_back(3,y1,x1,x2);
			state.emplace_back(1,y2,x1,x2);
		}
	}
	
	sort(state.begin(),state.end());
	
	recursive_set_2D_VG(0,s32(Y.size())-1,0,s32(state.size())-1,P2[lay],state,Y,1);
}
inline void build_2D_VG_point(const DataSet &data,std::vector<std::pair<u32,u32>> *P1,std::vector<std::pair<u32,u32>> *P2,std::vector<s64> &Px,std::vector<s64> &Py)
{
	std::vector<std::pair<u32,u32>> P3[11];
	std::vector<std::future<void>> taskX,taskY;
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		taskX.emplace_back(std::async(build_2D_VG_X_point,lay,ref(data),P1,P2,ref(Px),ref(Py)));
		taskY.emplace_back(std::async(build_2D_VG_Y_point,lay,ref(data),P1,P3,ref(Px),ref(Py)));
	}
	for(s32 i=0;i<data.metal_layers;++i)
	{
		taskX[i].wait();
		
		if(debug_mode)
		std::cerr<<"P2["<<i+1<<"].size(): "<<P2[i+1].size()<<endl;
	
		for(const auto &p:P2[i+1])
		{
			P1[i+1].emplace_back(p);
		}
		std::vector<std::pair<u32,u32>>().swap(P2[i+1]);
		
		taskY[i].wait();
		
		if(debug_mode)
		std::cerr<<"P3["<<i+1<<"].size(): "<<P3[i+1].size()<<endl;
	
		for(const auto &p:P3[i+1])
		{
			P1[i+1].emplace_back(p);
		}		
		std::vector<std::pair<u32,u32>>().swap(P3[i+1]);
		
		set_dis(P1[i+1]);
	}
	
	if(debug_mode)
	std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
}

void VisingGraph::build_beta(const DataSet &data,bool is_not_connect)
{	
	static const int LIMIT_LAYER = DataSet::LIMIT_LAYER;
	
	showclock("start");
	
	std::vector<std::pair<u32,u32>> P1[LIMIT_LAYER];
	std::vector<std::pair<u32,u32>> P2[LIMIT_LAYER];
	std::vector<std::pair<u32,u32>> P3[LIMIT_LAYER];
	std::vector<std::pair<u32,u32>> P4[LIMIT_LAYER];
	
	get_original_P1(P1,data.metal_layers,xLine,yLine,Px,Py,data);
	showclock("get_original_P1");
	
	//*
	if(debug_mode)
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::cerr<<"P1["<<lay<<"].size(): "<<P1[lay].size()<<endl;
	}
	std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	//*/
	
	project_point_on_all_layer_beta(1,data.metal_layers,data,P1,P2,Px,Py);// insert more point
	showclock("project_point_on_all_layer_beta");
	
	//*
	if(debug_mode)
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::cerr<<"P1["<<lay<<"].size(): "<<P1[lay].size()<<endl;
	}
	if(debug_mode) std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	//*/
	
	point_project_to_XYLine(P1,P2,P3,P4,data,xLine,yLine,Px,Py);
	showclock("point_project_to_XYLine");
	
	//*
	if(debug_mode)
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::cerr<<"P1["<<lay<<"].size(): "<<P1[lay].size()<<endl;
	}
	if(debug_mode) std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	//*/
	
	recursive_set_3d_VG_point(1,data.metal_layers,data,P1,Px,Py);
	showclock("recursive_set_3d_VG_point");
	//*
	if(debug_mode)
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::cerr<<"P1["<<lay<<"].size(): "<<P1[lay].size()<<endl;
	}
	if(debug_mode) std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	//*/
	
	build_2D_VG_point(data,P1,P2,Px,Py);
	showclock("build_2D_VG_point");
	
	point_project_to_XYLine(P1,P2,P3,P4,data,xLine,yLine,Px,Py);
	showclock("point_project_to_XYLine");
	
	//*
	if(debug_mode)
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::cerr<<"P1["<<lay<<"].size(): "<<P1[lay].size()<<endl;
	}
	if(debug_mode) std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	//*/
	
	std::vector<size_t> Pv[LIMIT_LAYER];
	
	put_the_point_number(data.metal_layers,V_set,P1,Pv);
	showclock("put_the_point_number");
	
	DST.init(V_set.size());
	
	std::vector<bool> ori_is_P(V_set.size());
	
	merge_same_via_point(DST,data,V_set,Px,Py);
	showclock("merge_same_via_point");
	
	set_shape_point(ori_is_P,data.metal_layers,xLine,yLine,V_set,DST);
	showclock("set_shape_point");
	
	merge_same_shape_point(DST,xLine,yLine,Pv,V_set,Px,Py,data.metal_layers);
	showclock("merge_same_shape_point");
	
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::vector<Statemant_2D_VG>().swap(xLine[lay]);
		std::vector<Statemant_2D_VG>().swap(yLine[lay]);
	}
	
	N = shrink_point(shrink_from,is_pinv,ori_is_P,DST,V_set.size());
	showclock("shrink_point");
	
	std::cerr<<N<<' '<<V_set.size()<<endl;
	
	edge.clear();
	edge.reserve(V_set.size()*12);
	showclock("edge.reserve(V_set.size()*12)");
	
	build_via_edge(data,shrink_from,Pv,Px,Py,edge,V_set);
	showclock("build_via_edge");
	
	build_2D_edge(data,shrink_from,Pv,Px,Py,edge,V_set);
	showclock("build_2D_edge");
	
	std::cerr<<"edge size: "<<edge.size()<<'\n';
	
	set_edge_and_graph(data.viacost,N,G,edge,shrink_from,Px,Py,V_set);
	showclock("set_edge_and_graph");
}