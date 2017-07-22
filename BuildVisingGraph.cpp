#include"DataLoader.h"
#include"BuildVisingGraph.h"
#include"SwapLineP1.h"
#include"SwapLineStatement.h"
#include"BinaryIndexTree.h"
#include"DisjoinSet.h"

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

inline void get_Line_swap_line(std::vector<Statemant_2D_VG> &line,std::vector<statementP1> &state,swape_line_P1 &SL,u32 l,u32 r)
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
		s32 type=st.inout>0?3:1;
		for(const auto &seg:SL.segments)
		{
			line.emplace_back(type,st.a,seg.first,seg.second,st.seg_type);
		}
		SL.segments.clear();
	}
}

inline void get_xyLine_singal_layer(s32 lay,std::vector<Statemant_2D_VG> &xLine,std::vector<Statemant_2D_VG> &yLine,const DataSet &data,const std::vector<s64> &Px,const std::vector<s64> &Py)
{
	std::vector<statementP1> Xstate,Ystate;
	swape_line_P1 SL;
	
	for(auto o:data.Obstacles[lay])
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
	for(auto o:data.RoutedShape[lay])
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
	
	get_Line_swap_line(xLine,Xstate,SL,0,Py.size());
	get_Line_swap_line(yLine,Ystate,SL,0,Px.size());
}

inline void get_PxPy(std::vector<s64> &Px,std::vector<s64> &Py,const DataSet &data,std::vector<Statemant_2D_VG> *xLine,std::vector<Statemant_2D_VG> *yLine)
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
		
		if(x1<=x2&&y1<y2)
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
	if(P2==nullptr) return;
	
	for(const auto &p:S)
	{
		P2[la2].emplace_back(p);
	}
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

/*
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
}
*/
	
void one_way_point_project_swap_line(const Statemant_2D_VG &st,std::set<u32> &ST,std::vector<std::pair<u32,u32>> &S2,u32 L,u32 R,bool is_rev=0)
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
		
		while(it_l!=it_r)
		{
			auto tmp=it_l++;
			if(!is_rev)
			{
				S2.emplace_back(st.a,*tmp);
			}
			else
			{
				S2.emplace_back(*tmp,st.a);
			}
			ST.erase(tmp);
		}
	}
}

void single_layer_point_project(const std::vector<std::pair<u32,u32>> &S,std::vector<Statemant_2D_VG> &xLine,std::vector<Statemant_2D_VG> &yLine,std::vector<std::pair<u32,u32>> &S2,u32 x1,u32 x2,u32 y1,u32 y2)
{
	std::vector<Statemant_2D_VG> Xstate=xLine;
	std::vector<Statemant_2D_VG> Ystate=yLine;
	
	for(const auto &p:S)
	{
		Xstate.emplace_back(2,p.first,p.second);
		Ystate.emplace_back(2,p.second,p.first);
	}
	
	sort(Xstate.begin(),Xstate.end());
	sort(Ystate.begin(),Ystate.end());
	
	std::set<u32> ST,RST;
	
	for(size_t i=0;i<Xstate.size();++i)
	{
		size_t ri=Xstate.size()-i-1;
		//if(Xstate[i].seg_type!='B'||Xstate[i].type!=3)
			one_way_point_project_swap_line(Xstate[i],ST,S2,x1,x2);
		//if(Xstate[ri].seg_type!='B'||Xstate[ri].type!=1)
			one_way_point_project_swap_line(Xstate[ri],RST,S2,x1,x2);
	}
	
	ST.clear();
	RST.clear();
	
	for(size_t i=0;i<Ystate.size();++i)
	{
		size_t ri=Ystate.size()-i-1;
		//if(Ystate[i].seg_type!='B'||Ystate[i].type!=3)
			one_way_point_project_swap_line(Ystate[i],ST,S2,y1,y2,1);
		//if(Ystate[ri].seg_type!='B'||Ystate[ri].type!=1)
			one_way_point_project_swap_line(Ystate[ri],RST,S2,y1,y2,1);
	}
	
}

inline void point_project_to_XYLine(std::vector<std::pair<u32,u32>> *P1,std::vector<std::pair<u32,u32>> *P2,const DataSet &data,std::vector<Statemant_2D_VG> *xLine,std::vector<Statemant_2D_VG> *yLine,const std::vector<s64> &Px,const std::vector<s64> &Py)
{
	u32 x1=get_dis(Px,data.boundary.first.x);
	u32 x2=get_dis(Px,data.boundary.second.x);
	u32 y1=get_dis(Py,data.boundary.first.y);
	u32 y2=get_dis(Py,data.boundary.second.y);
	
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		/*
		xLine[lay].emplace_back(1,x1,y1,y2,'B');
		xLine[lay].emplace_back(3,x2,y1,y2,'B');
		yLine[lay].emplace_back(1,y1,x1,x2,'B');
		yLine[lay].emplace_back(3,y2,x1,x2,'B');
		*/
		
		single_layer_point_project(P1[lay],xLine[lay],yLine[lay],P2[lay],x1,x2,y1,y2);
		for(const auto &p:P1[lay])
		{
			P2[lay].emplace_back(p);
		}
		set_dis(P2[lay]);
		P1[lay].clear();
		for(const auto &p:P2[lay])
		{
			if(x1<=p.first&&p.first<=x2&&y1<=p.second&&p.second<=y2)
				P1[lay].emplace_back(p);
		}
		P2[lay].clear();
	}
	
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
			ST.emplace(state[sl_mid].b1);
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
			ST.emplace(state[sr_mid].b1);
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

inline void build_2D_VG_X_point(s32 lay,const DataSet &data,std::vector<std::pair<u32,u32>> *P1,std::vector<std::pair<u32,u32>> *P2,std::vector<s64> &Px,std::vector<s64> &Py)
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
		s32 y2=s32(get_upper_dis(Py,o.second.y))-1;
		
		if(x1<=x2&&y1<=y2)
		{
			state.emplace_back(3,x1,y1,y2);
			state.emplace_back(1,x2,y1,y2);
		}
	}
	
	sort(state.begin(),state.end());
	
	recursive_set_2D_VG(0,s32(X.size())-1,0,s32(state.size())-1,P2[lay],state,X);
}

inline void build_2D_VG_Y_point(s32 lay,const DataSet &data,std::vector<std::pair<u32,u32>> *P1,std::vector<std::pair<u32,u32>> *P2,std::vector<s64> &Px,std::vector<s64> &Py)
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
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		build_2D_VG_X_point(lay,data,P1,P2,Px,Py);
		
		std::cerr<<"P2["<<lay<<"].size(): "<<P2[lay].size()<<endl;
		
		build_2D_VG_Y_point(lay,data,P1,P2,Px,Py);
		
		std::cerr<<"  P2["<<lay<<"].size(): "<<P2[lay].size()<<endl;
		
		for(const auto &p:P2[lay])
		{
			P1[lay].emplace_back(p);
		}
		set_dis(P1[lay]);
		
		P2[lay]=std::vector<std::pair<u32,u32>>();
	}
	std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
}

inline void put_the_point_number(s32 metal_layers,std::vector<point3D> &V_set,std::vector<std::tuple<u32,u32,size_t>> *V,std::vector<std::pair<u32,u32>> *P1)
{
	for(s32 lay=1;lay<=metal_layers;++lay)
	{
		for(const auto &p:P1[lay])
		{
			V_set.emplace_back(p.first,p.second,lay);
		}
	}
	set_dis(V_set);
	
	for(s32 lay=1;lay<=metal_layers;++lay)
	{
		for(const auto &p:P1[lay])
		{
			V[lay].emplace_back(p.first,p.second,get_dis(V_set,point3D(p.first,p.second,lay)));
		}
		P1[lay]=std::vector<std::pair<u32,u32>>();
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
			
			if(!(V_set[p1]==P3D1))std::cerr<<"Error V_set[p1]!=P3D1\n";
			if(!(V_set[p2]==P3D2))std::cerr<<"Error V_set[p2]!=P3D2\n";
			
			DST.U(p1,p2);
		}
	}
}

inline void set_shape_point(std::vector<bool> &ori_is_P,s32 metal_layers,std::vector<Statemant_2D_VG> *xLine,std::vector<Statemant_2D_VG> *yLine,const std::vector<point3D> &V_set,DisjoinSet &DST)
{
	for(s32 lay=1;lay<=metal_layers;++lay)
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
			
			if(p1<V_set.size()&&V_set[p1]==P3D1)
				ori_is_P[p1]=true;
			else std::cerr<<"Error V_set[p1]!=P3D1 .. X\n";
			if(p2<V_set.size()&&V_set[p2]==P3D2)
				ori_is_P[p2]=true;
			else std::cerr<<"Error V_set[p2]!=P3D2 .. X\n";
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
			
			if(p1<V_set.size()&&V_set[p1]==P3D1)
				ori_is_P[p1]=true;
			else cout<<"Error V_set[p1]!=P3D1 .. Y\n";
			if(p2<V_set.size()&&V_set[p2]==P3D2)
				ori_is_P[p2]=true;
			else cout<<"Error V_set[p2]!=P3D2 .. Y\n";
		}
	}
}

inline void merge_shape_point_swap_line(s32 lay,DisjoinSet &DST,std::vector<Statemant_2D_VG> &state,const std::vector<point3D> &V_set,u32 bit_size,bool is_rev=0)
{
	sort(state.begin(),state.end());
	
	std::map<u32,u32> ST;
	
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
					point3D P3D1(it->second,it->first,lay);
					point3D P3D2(st.a,st.b1,lay);
					
					if(is_rev)
					{
						std::swap(P3D1.x,P3D1.y);
						std::swap(P3D2.x,P3D2.y);
					}
					
					size_t p1=get_dis(V_set,P3D1);
					size_t p2=get_dis(V_set,P3D2);
					
					if(!(V_set[p1]==P3D1))std::cerr<<"Error V_set[p1]!=P3D1\n";
					if(!(V_set[p2]==P3D2))std::cerr<<"Error V_set[p2]!=P3D2\n";
					
					DST.U(p1,p2);
				}
				ST[st.b1]=st.a;
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
			while(it_l!=it_r)
			{
				auto tmp=it_l++;
				ST.erase(tmp);
			}
		}
	}
	
	ST.clear();
	BIST.init(bit_size+2);
	
	for(s32 jj=s32(state.size())-1;jj>=0;--jj)
	{
		auto &st=state[jj];
		if(st.type==2)
		{
			if(BIST.get_sum(st.b1+1)||BIST.get_sum(st.b1))
			{
				auto it=ST.find(st.b1);
				if(it!=ST.end())
				{
					point3D P3D1(it->second,it->first,lay);
					point3D P3D2(st.a,st.b1,lay);
					
					if(is_rev)
					{
						std::swap(P3D1.x,P3D1.y);
						std::swap(P3D2.x,P3D2.y);
					}
					
					size_t p1=get_dis(V_set,P3D1);
					size_t p2=get_dis(V_set,P3D2);
					
					if(!(V_set[p1]==P3D1))std::cerr<<"Error V_set[p1]!=P3D1\n";
					if(!(V_set[p2]==P3D2))std::cerr<<"Error V_set[p2]!=P3D2\n";
					
					DST.U(p1,p2);
				}
				ST[st.b1]=st.a;
			}
		}
		else
		{
			if(st.type==1)
			{
				BIST.add(st.b1+1,-1);
				BIST.add(st.b2+1,1);
			}
			else
			{
				BIST.add(st.b1+1,1);
				BIST.add(st.b2+1,-1);
			}
			if(st.type==3)continue;
			auto it_l=ST.lower_bound(st.b1);
			auto it_r=ST.upper_bound(st.b2);
			while(it_l!=it_r)
			{
				auto tmp=it_l++;
				ST.erase(tmp);
			}
		}
	}
}

inline void merge_same_shape_point(DisjoinSet &DST,std::vector<Statemant_2D_VG> *xLine,std::vector<Statemant_2D_VG> *yLine,std::vector<std::tuple<u32,u32,size_t>> *V,const std::vector<point3D> &V_set,const std::vector<s64> &Px,const std::vector<s64> &Py,s32 metal_layers)
{
	for(s32 lay=1;lay<=metal_layers;++lay)
	{
		std::vector<Statemant_2D_VG> Xstate;
		std::vector<Statemant_2D_VG> Ystate;
				
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
		
		for(const auto &p:V[lay])
		{
			Xstate.emplace_back(2,std::get<0>(p),std::get<1>(p));
			Ystate.emplace_back(2,std::get<1>(p),std::get<0>(p));
		}
		
		merge_shape_point_swap_line(lay,DST,Xstate,V_set,Py.size());
		merge_shape_point_swap_line(lay,DST,Ystate,V_set,Px.size(),1);
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
	
	is_pinv.resize(tmp.size());
	
	for(size_t i=0;i<sz;++i)
	{
		shrink_from[i]=get_dis(tmp,shrink_from[i]);
		is_pinv[shrink_from[i]] = is_pinv[shrink_from[i]] || ori_is_P[i];
	}
	
	return tmp.size();
}

inline void build_2D_edge_swape_line(s32 lay,const std::vector<size_t> &shrink_from,std::vector<Statemant_2D_VG> &state,const std::vector<point3D> &V_set,std::vector<Edge> &edge,bool is_rev=0)
{
	sort(state.begin(),state.end());
	
	std::map<u32,u32> ST;
	
	for(const auto &st:state)
	{
		if(st.type==2)
		{
			auto it=ST.find(st.b1);
			if(it!=ST.end())
			{
				point3D P3D1(it->second,it->first,lay);
				point3D P3D2(st.a,st.b1,lay);
				
				u8 edge_type='H';
				if(is_rev)
				{
					std::swap(P3D1.x,P3D1.y);
					std::swap(P3D2.x,P3D2.y);
					edge_type='V';
				}
				
				size_t p1=get_dis(V_set,P3D1);
				size_t p2=get_dis(V_set,P3D2);
				
				if(!(V_set[p1]==P3D1))std::cerr<<"Error V_set[p1]!=P3D1\n";
				if(!(V_set[p2]==P3D2))std::cerr<<"Error V_set[p2]!=P3D2\n";
				
				if(shrink_from[p1]!=shrink_from[p2])
				{
					edge.emplace_back(p1,p2,edge_type);
					edge.emplace_back(p2,p1,edge_type);
				}
			}
			ST[st.b1]=st.a;
		}
		else
		{
			auto it_l=ST.lower_bound(st.b1);
			auto it_r=ST.upper_bound(st.b2);
			while(it_l!=it_r)
			{
				auto tmp=it_l++;
				ST.erase(tmp);
			}
		}
	}
	
}

inline void build_2D_edge(const DataSet &data,const std::vector<size_t> &shrink_from,std::vector<std::tuple<u32,u32,size_t>> *V,std::vector<s64> &Px,std::vector<s64> &Py,std::vector<Edge> &edge,const std::vector<point3D> &V_set)
{
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::vector<Statemant_2D_VG> Xstate;
		std::vector<Statemant_2D_VG> Ystate;
		
		for(const auto &p:V[lay])
		{
			Xstate.emplace_back(2,std::get<0>(p),std::get<1>(p));
			Ystate.emplace_back(2,std::get<1>(p),std::get<0>(p));
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
		
		build_2D_edge_swape_line(lay,shrink_from,Xstate,V_set,edge);
		build_2D_edge_swape_line(lay,shrink_from,Ystate,V_set,edge,1);
	}
}

inline void build_via_edge_swap_line(s32 la1,std::vector<Statemant_2D_VG> &state_la2,std::vector<std::tuple<u32,u32,size_t>> *V,std::vector<std::tuple<u32,u32,s32>> &P,std::map<std::pair<u32,u32>,s32> &res,u32 bit_size)
{
	for(const auto &p:P)
	{
		state_la2.emplace_back(2,std::get<0>(p),std::get<1>(p),std::get<2>(p));
	}
	
	P.clear();
	res.clear();
	
	for(const auto &p:V[la1])
	{
		state_la2.emplace_back(2,std::get<0>(p),std::get<1>(p),la1);
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
					res[{st.a,st.b1}]=st.b2;
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

inline void build_via_edge(const DataSet &data,const std::vector<size_t> &shrink_from,std::vector<std::tuple<u32,u32,size_t>> *V,std::vector<s64> &Px,std::vector<s64> &Py,std::vector<Edge> &edge,const std::vector<point3D> &V_set)
{
	std::vector<std::tuple<u32,u32,s32>> P;
	std::map<std::pair<u32,u32>,s32> res;
	
	for(s32 lay=2;lay<=data.metal_layers;++lay)
	{
		std::vector<Statemant_2D_VG> state;
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
		
		build_via_edge_swap_line(lay-1,state,V,P,res,Py.size());
				
		for(const auto &p:V[lay])
		{
			auto it=res.find({std::get<0>(p),std::get<1>(p)});
			if(it!=res.end())
			{
				point3D P3D1(it->first.first,it->first.second,it->second);
				size_t p1=get_dis(V_set,P3D1);
				size_t p2=std::get<2>(p);
				if(shrink_from[p1]!=shrink_from[p2])
				{
					edge.emplace_back(p1,p2,'Z');
					edge.emplace_back(p2,p1,'Z');
				}
				res.erase(it);
			}
		}
		
		for(const auto &p:res)
		{
			P.emplace_back(p.first.first,p.first.second,p.second);
		}
		
	}
}

inline void set_edge_and_graph(s64 viacost,size_t N,std::vector<std::vector<size_t>> &G,std::vector<Edge> &edge,const std::vector<size_t> &shrink_from,std::vector<s64> &Px,std::vector<s64> &Py,const std::vector<point3D> &V_set)
{
	G.clear();
	G.resize(N);
	
	for(size_t i=0;i<edge.size();++i)
	{
		auto &e=edge[i];
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
		
		e.u=shrink_from[e.ori_u];
		e.v=shrink_from[e.ori_v];
		
		if(e.u==e.v)std::cerr<<"Error!\n";
		
		G[e.u].emplace_back(i);
		
	}
}

/*
void find_shape_point_swap_line(u32 bit_size,std::vector<Statemant_2D_VG> &state,std::vector<point3D> &P,s32 lay)
{
	sort(state.begin(),state.end());
		
	BIT BIST;
	BIST.init(bit_size+2);
	
	for(const auto &st:state)
	{
		if(st.type==2)
		{
			if(BIST.get_sum(st.b1+1)||BIST.get_sum(st.b1))
			{
				P.emplace_back(st.a,st.b1,lay);
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
		}
	}
}

void find_shape_point(std::vector<point3D> &P,std::vector<Statemant_2D_VG> *xLine,std::vector<Statemant_2D_VG> *yLine,std::vector<std::tuple<u32,u32,size_t>> *V,const std::vector<point3D> &V_set,const std::vector<s64> &Px,const std::vector<s64> &Py,s32 metal_layers)
{
	for(s32 lay=1;lay<=metal_layers;++lay)
	{
		std::vector<Statemant_2D_VG> Xstate;
		std::vector<Statemant_2D_VG> Ystate;
				
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
		
		for(const auto &p:V[lay])
		{
			Xstate.emplace_back(2,std::get<0>(p),std::get<1>(p));
			Ystate.emplace_back(2,std::get<1>(p),std::get<0>(p));
		}
		find_shape_point_swap_line(Py.size(),Xstate,P,lay);
	}
}
*/

void VisingGraph::build(const DataSet &data)
{	
	static const int LIMIT_LAYER = DataSet::LIMIT_LAYER;
	
	std::vector<Statemant_2D_VG> xLine[LIMIT_LAYER];
	std::vector<Statemant_2D_VG> yLine[LIMIT_LAYER];
	
	get_original_PxPy(Px,Py,data);
	
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		get_xyLine_singal_layer(lay,xLine[lay],yLine[lay],data,Px,Py);
	}
	/*
	for(s32 lay=1;lay<=data.metal_layers;++lay){
		for(auto l:yLine[lay])
		{
			if(l.seg_type=='S'){
				cout<<"H-line M"<<lay<<" ("<<Px[l.b1]<<","<<Py[l.a]<<") ("<<Px[l.b2]<<","<<Py[l.a]<<")\n";
			}
		}
	}
	//*/
	get_PxPy(Px,Py,data,xLine,yLine);
	
	std::vector<std::pair<u32,u32>> P1[LIMIT_LAYER];
	std::vector<std::pair<u32,u32>> P2[LIMIT_LAYER];
	
	get_original_P1(P1,data.metal_layers,xLine,yLine,Px,Py,data);
	
	//*
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::cerr<<"P1["<<lay<<"].size(): "<<P1[lay].size()<<endl;
	}
	std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	//*/
	
	recursive_set_3d_VG_point(1,data.metal_layers,data,P1,Px,Py);
	//*
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::cerr<<"P1["<<lay<<"].size(): "<<P1[lay].size()<<endl;
	}
	std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	//*/
	
	point_project_to_XYLine(P1,P2,data,xLine,yLine,Px,Py);
	
	/*
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		for(const auto &p:P1[lay])
		{
			if(lay==1)
			cout<<"Via V"<<lay<<" ("<<Px[p.first]<<","<<Py[p.second]<<")\n";
		}
	}
	//*/
	
	//*
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::cerr<<"P1["<<lay<<"].size(): "<<P1[lay].size()<<endl;
	}
	std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	//*/
	
	recursive_set_3d_VG_point(1,data.metal_layers,data,P1,Px,Py);
	//*
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::cerr<<"P1["<<lay<<"].size(): "<<P1[lay].size()<<endl;
	}
	std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	//*/
	
	build_2D_VG_point(data,P1,P2,Px,Py);
	
	point_project_to_XYLine(P1,P2,data,xLine,yLine,Px,Py);//add more point, delete OK
	
	/*
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		for(const auto &p:P1[lay])
		{
			if(lay==1)
			cout<<"Via V"<<lay<<" ("<<Px[p.first]<<","<<Py[p.second]<<")\n";
		}
	}
	//*/
	
	//*
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::cerr<<"P1["<<lay<<"].size(): "<<P1[lay].size()<<endl;
	}
	std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	//*/
	
	std::vector<std::tuple<u32,u32,size_t>> V[LIMIT_LAYER];
	
	put_the_point_number(data.metal_layers,V_set,V,P1);
	
	//*
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		std::cerr<<"V["<<lay<<"].size(): "<<V[lay].size()<<endl;
	}
	std::cerr<<"wwwwwwwwwwwwwwwwwwwwwwwwwwww\n";
	//*/
	
	DisjoinSet DST(V_set.size());
	
	std::vector<bool> ori_is_P(V_set.size());
	
	merge_same_via_point(DST,data,V_set,Px,Py);
	set_shape_point(ori_is_P,data.metal_layers,xLine,yLine,V_set,DST);
	merge_same_shape_point(DST,xLine,yLine,V,V_set,Px,Py,data.metal_layers);
	
	/*
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		for(const auto &p:V[lay])
		{
			if(DST.size(std::get<2>(p))>1&&lay==1)
			cout<<"Via V"<<lay<<" ("<<Px[std::get<0>(p)]<<","<<Py[std::get<1>(p)]<<")\n";
		}
	}
	//*/
	
	N = shrink_point(shrink_from,is_pinv,ori_is_P,DST,V_set.size());
	std::cerr<<N<<' '<<V_set.size()<<endl;
	/*
	19720
	*/
	
	/*
	size_t gggg;
	std::cin>>gggg;
	for(s32 lay=1;lay<=data.metal_layers;++lay)
	{
		for(const auto &p:V[lay])
		{
			if(shrink_from[std::get<2>(p)]==gggg)
			cout<<"Via V"<<lay<<" ("<<Px[std::get<0>(p)]<<","<<Py[std::get<1>(p)]<<")\n";
		}
	}
	//*/
	
	build_2D_edge(data,shrink_from,V,Px,Py,edge,V_set);
	
	/*
	for(size_t i=0;i<edge.size();i+=2)
	{
		cout<<edge[i].type<<"-line M"<<V_set[edge[i].ori_u].layer<<" ("<<Px[V_set[edge[i].ori_u].x]<<","<<Py[V_set[edge[i].ori_u].y]<<") ("<<Px[V_set[edge[i].ori_v].x]<<","<<Py[V_set[edge[i].ori_v].y]<<")\n";
	}
	//*/
	
	build_via_edge(data,shrink_from,V,Px,Py,edge,V_set);
	
	/*
	for(size_t i=0;i<edge.size();i+=2)
	{
		if(edge[i].type=='Z'){
			for(auto lay=V_set[edge[i].ori_u].layer;lay<V_set[edge[i].ori_v].layer;++lay)
			cout<<"Via V"<<lay<<" ("<<Px[V_set[edge[i].ori_u].x]<<","<<Py[V_set[edge[i].ori_u].y]<<")\n";
		}
		else{
			cout<<edge[i].type<<"-line M"<<V_set[edge[i].ori_u].layer<<" ("<<Px[V_set[edge[i].ori_u].x]<<","<<Py[V_set[edge[i].ori_u].y]<<") ("<<Px[V_set[edge[i].ori_v].x]<<","<<Py[V_set[edge[i].ori_v].y]<<")\n";
		}
	}
	//*/
	
	
	set_edge_and_graph(data.viacost,N,G,edge,shrink_from,Px,Py,V_set);
	//*
	s64 cost=0;
	for(size_t i=0;i<edge.size();i+=2)
	{
		cost+=edge[i].cost;
	}
	std::cerr<<cost<<endl;
	
	u32 cnt=0;
	std::vector<size_t> tmd,tmd2;
	for(size_t i=0;i<is_pinv.size();++i)if(is_pinv[i])tmd.push_back(i),++cnt;
	/*
	for(auto p:V[1]){
		if(get_dis(tmd,std::get<2>(p))<tmd.size()&&tmd[get_dis(tmd,std::get<2>(p))]==std::get<2>(p))tmd2.push_back(std::get<2>(p));
	}
	set_dis(tmd2);
	for(auto i:tmd2)std::cerr<<i<<endl;
	//*/
	std::cerr<<cnt<<endl;
	//*/
}