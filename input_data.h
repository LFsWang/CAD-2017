#include<cstdio>
#include<cstring>
#include<vector>
#include<algorithm>
#include<sstream>

using std::vector;
using std::stringstream;
using std::max;

struct input_data{
	struct __RoutedShape{
		unsigned Layer;
		long long LLx,LLy,URx,URy;
		__RoutedShape(unsigned Layer,long long LLx,long long LLy,long long URx,long long URy):
			Layer(Layer),LLx(LLx),LLy(LLy),URx(URx),URy(URy){}
	};
	struct __RoutedVia{
		unsigned Layer;
		long long x,y;
		__RoutedVia(unsigned Layer,long long x,long long y):
			Layer(Layer),x(x),y(y){}
	};
	struct __Obstacle{
		unsigned Layer;
		long long LLx,LLy,URx,URy;
		__Obstacle(unsigned Layer,long long LLx,long long LLy,long long URx,long long URy):
			Layer(Layer),LLx(LLx),LLy(LLy),URx(URx),URy(URy){}
	};
	
	vector<long long> dis;
	int ViaCost, Spacing;
	long long LLx,LLy,URx,URy;
	unsigned MetalLayers;
	vector<__RoutedShape> RoutedShape;
	vector<__RoutedVia> RoutedVia;
	vector<__Obstacle> Obstacle;
	
	void init_by_stdin();
	void print()const;
	unsigned get_dis(long long d)const;
};
void input_data::init_by_stdin(){
	dis = {LLONG_MAX, LLONG_MIN};
	
	RoutedShape.clear();
	RoutedVia.clear();
	Obstacle.clear();
	
	scanf("%*s %*s %d",&ViaCost);
	dis.emplace_back(ViaCost);
	scanf("%*s %*s %d",&Spacing);
	scanf("%*s %*s (%lld,%lld) (%lld,%lld)",&LLx,&LLy,&URx,&URy);
	
	dis.emplace_back(LLx);
	dis.emplace_back(LLy);
	dis.emplace_back(URx);
	dis.emplace_back(URy);
	
	scanf("%*s %*s %u",&MetalLayers);
	
	unsigned RoutedShapes;
	unsigned RoutedVias;
	unsigned Obstacles;
	scanf("%*s %*s %u",&RoutedShapes);
	scanf("%*s %*s %u",&RoutedVias);
	scanf("%*s %*s %u",&Obstacles);
	
	int n=RoutedShapes+RoutedVias+Obstacles;
	
	for(char type[15],M[5];n--;){
		scanf("%s%s",type,M);
		unsigned Layer;
		long long lx,ly,ux,uy;
		stringstream ss(M+1);
		ss>>Layer;
		if(!strcmp(type,"RoutedShape")){
			scanf(" (%lld,%lld) (%lld,%lld)",&lx,&ly,&ux,&uy);
			dis.emplace_back(lx);
			dis.emplace_back(ly);
			dis.emplace_back(ux);
			dis.emplace_back(uy);
			RoutedShape.emplace_back(Layer,lx,ly,ux,uy);
		}else if(!strcmp(type,"RoutedVia")){
			scanf(" (%lld,%lld)",&lx,&ly);
			dis.emplace_back(lx);
			dis.emplace_back(ly);
			RoutedVia.emplace_back(Layer,lx,ly);
		}else{
			scanf(" (%lld,%lld) (%lld,%lld)",&lx,&ly,&ux,&uy);
			Obstacle.emplace_back(Layer,lx-Spacing,ly-Spacing,ux+Spacing,uy+Spacing);
			dis.emplace_back(lx-Spacing);
			dis.emplace_back(ly-Spacing);
			dis.emplace_back(ux+Spacing);
			dis.emplace_back(uy+Spacing);
		}
	}
	sort(dis.begin(),dis.end());
	dis.resize(unique(dis.begin(),dis.end())-dis.begin());
}
void input_data::print()const{
	long long d = 100;
	
	printf("ViaCost = %d\n",ViaCost);
	printf("Spacing = 0\n");
	printf("Boundary = (%lld,%lld) (%lld,%lld)\n",LLx,LLy,URx+d,URy+d);
	printf("#MetalLayers = %u\n",MetalLayers);
	printf("#RoutedShapes = %u\n",unsigned(RoutedShape.size()));
	printf("#RoutedVias = %u\n",unsigned(RoutedVia.size()));
	printf("#Obstacles = %u\n",unsigned(Obstacle.size()));
	
	for(const auto &i:RoutedShape){
		printf("RoutedShape M%u (%lld,%lld) (%lld,%lld)\n",i.Layer,i.LLx,i.LLy,i.URx,i.URy);
	}
	for(const auto &i:RoutedVia){
		printf("RoutedVia M%u (%lld,%lld)\n",i.Layer,i.x,i.y);
	}
	for(const auto &i:Obstacle){
		printf("Obstacle M%u (%lld,%lld) (%lld,%lld)\n",i.Layer,max(i.LLx,LLx),max(i.LLy,LLy),max(i.URx,LLx),max(i.URy,LLy));
	}
}
unsigned input_data::get_dis(long long d)const{
	return lower_bound(dis.begin(),dis.end(),d) - dis.begin();
}
