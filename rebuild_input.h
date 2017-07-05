#include"input_data.h"
#include"swape_line.h"
namespace rebuild_input_by_swape_line{
	struct statement{
		unsigned x,y1,y2;//y2 > y1
		char type;// S or O
		int inout;//-1 or 1
		statement(unsigned x,unsigned y1,unsigned y2,char type,int inout):
				x(x),y1(y1),y2(y2),type(type),inout(inout){}
		bool operator<(const statement &b)const{
			return make_tuple(x,inout) < make_tuple(b.x,b.inout) ;
		}
	};
	
	vector<vector<statement>> get_statement(const input_data &input){
		vector<vector<statement>> res(input.MetalLayers+1);
		for(const auto &s:input.RoutedShape){
			unsigned lx=input.get_dis(s.LLx);
			unsigned ly=input.get_dis(s.LLy);
			unsigned ux=input.get_dis(s.URx);
			unsigned uy=input.get_dis(s.URy);
			res[s.Layer].emplace_back(lx,ly,uy,'S',1);
			res[s.Layer].emplace_back(ux,ly,uy,'S',-1);
		}
		for(const auto &o:input.Obstacle){
			unsigned lx=input.get_dis(o.LLx);
			unsigned ly=input.get_dis(o.LLy);
			unsigned ux=input.get_dis(o.URx);
			unsigned uy=input.get_dis(o.URy);
			res[o.Layer].emplace_back(lx,ly,uy,'O',1);
			res[o.Layer].emplace_back(ux,ly,uy,'O',-1);
		}
		for(auto &r:res)sort(r.begin(),r.end());
		return res;
	}
	input_data rebuild_input(const input_data &input){
		vector<vector<statement>> res = get_statement(input);
		input_data output;
		swape_line SL;
		unsigned n = input.dis.size();
		for(size_t layer = 1;layer<res.size();++layer){
			SL.init(n);
			for(const auto &st:res[layer]){
				auto LT=SL.findL(st.y1,0,0,n-1,1);
				auto RT=SL.findR(st.y2,0,0,n-1,1);
				
				SL.segments.clear();
				
				if( st.inout > 0 ){
					if(!SL.is_same(st.y1,st.y2,st.type,0,n-1,1)){
						SL.segments.emplace_back(get<2>(LT),st.y1,get<1>(LT),get<0>(LT));
						SL.find_seg(st.y1,st.y2,0,n-1,1);
						SL.segments.emplace_back(st.y2,get<2>(RT),get<1>(RT),get<0>(RT));
						SL.set_deep(get<0>(SL.segments.front()),get<1>(SL.segments.back()),st.x,0,n-1,1);
					}//else puts("GG");
					SL.insert(st.y1,st.y2,st.inout,st.type,0,n-1,1);
				}else{
					SL.segments.emplace_back(get<2>(LT),st.y1,get<1>(LT),get<0>(LT));
					SL.find_seg(st.y1,st.y2,0,n-1,1);
					SL.segments.emplace_back(st.y2,get<2>(RT),get<1>(RT),get<0>(RT));
					SL.insert(st.y1,st.y2,st.inout,st.type,0,n-1,1);
					if(!SL.is_same(st.y1,st.y2,st.type,0,n-1,1)){
						SL.set_deep(get<0>(SL.segments.front()),get<1>(SL.segments.back()),st.x,0,n-1,1);
					}else SL.segments.clear();//,puts("GG");
				}
				
				vector<tuple<unsigned,unsigned,char,int>> segments;
				//(l,r) type deep
				
				for(auto &i:SL.segments){
					//printf("	(%u,%u) %c %d\n",get<0>(i),get<1>(i),get<2>(i),get<3>(i));
					if(segments.size()&&get<2>(segments.back())==get<2>(i)&&get<3>(segments.back())==get<3>(i)&&get<1>(segments.back())==get<0>(i)){
						get<1>(segments.back())=get<1>(i);
					}else segments.emplace_back(i);
				}
				for(auto &i:segments){
					if(st.x==unsigned(get<3>(i)))continue;
					long long lx=get<3>(i),ly=get<0>(i),ux=st.x,uy=get<1>(i);
					lx=input.dis[lx];
					ly=input.dis[ly];
					ux=input.dis[ux];
					uy=input.dis[uy];
					switch(get<2>(i)){
						case 'S':{
							output.RoutedShape.emplace_back(layer,lx,ly,ux,uy);
							break;
						}
						case 'O':{
							output.Obstacle.emplace_back(layer,lx,ly,ux,uy);
							break;
						}
						default:{
							
						}
					}
				}
			}
		}
		output.dis=input.dis;
		output.LLx=input.LLx;
		output.LLy=input.LLy;
		output.URx=input.URx;
		output.URy=input.URy;
		output.ViaCost=input.ViaCost;
		output.Spacing=input.Spacing;
		output.MetalLayers=input.MetalLayers;
		output.RoutedVia=input.RoutedVia;
		return output;
	}
}
