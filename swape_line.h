#include<vector>
#include<utility>
#include<algorithm>

using std::vector;
using std::tuple;
using std::make_tuple;
using std::get;
using std::max;
using std::min;

struct swape_line{
	struct node{
		unsigned obstacle_tag,shape_tag;
		unsigned obstacle_sum,shape_sum;
		int deep_tag,deep;
		node():obstacle_tag(0),shape_tag(0),obstacle_sum(0),shape_sum(0),deep_tag(-1),deep(0){};
	};
	
	unsigned n;
	vector<node> st;
	vector<tuple<unsigned,unsigned,char,int>> segments;
	//(l,r) type deep
	
	
	unsigned get_obstacle_sum(unsigned l,unsigned r,unsigned d);
	unsigned get_shape_sum(unsigned l,unsigned r,unsigned d);
	int get_deep(unsigned d);
	void deep_up(unsigned d);
	void up(unsigned l,unsigned r,unsigned d);
	void down(unsigned d);
	
	void init(unsigned _n);
	
	void insert(unsigned a,unsigned b,int w,char type,unsigned l,unsigned r,unsigned d);
	void set_deep(unsigned a,unsigned b,int dep,unsigned l,unsigned r,unsigned d);
	tuple<int,char,unsigned> findL(unsigned x,char type,unsigned l,unsigned r,unsigned d);
	tuple<int,char,unsigned> findR(unsigned x,char type,unsigned l,unsigned r,unsigned d);
	bool is_same(unsigned a,unsigned b,char type,unsigned l,unsigned r,unsigned d);
	void find_seg(unsigned a,unsigned b,unsigned l,unsigned r,unsigned d);
};

unsigned swape_line::get_obstacle_sum(unsigned l,unsigned r,unsigned d){
	if(st[d].obstacle_tag) return r - l;
	return st[d].obstacle_sum;
}
unsigned swape_line::get_shape_sum(unsigned l,unsigned r,unsigned d){
	if(st[d].shape_tag) return r - l;
	return st[d].shape_sum;
}
int swape_line::get_deep(unsigned d){
	return st[d].deep_tag != -1 ? st[d].deep_tag : st[d].deep;
}
void swape_line::deep_up(unsigned d){
	if(st[d].deep_tag == -1) st[d].deep = get_deep(d*2);
}
void swape_line::up(unsigned l,unsigned r,unsigned d){
	unsigned mid = (l+r)/2;
	st[d].obstacle_sum = get_obstacle_sum(l,mid,d*2) + get_obstacle_sum(mid,r,d*2+1);
	st[d].shape_sum = get_shape_sum(l,mid,d*2) + get_shape_sum(mid,r,d*2+1);
	deep_up(d);
}
void swape_line::down(unsigned d){
	if(~st[d].deep_tag){
		st[d*2].deep_tag = st[d].deep_tag;
		st[d*2+1].deep_tag = st[d].deep_tag;
		st[d].deep_tag = -1;
	}
}
void swape_line::init(unsigned _n){
	st=vector<node>( (n = _n) * 4 );
	segments.clear();
}
void swape_line::insert(unsigned a,unsigned b,int w,char type,unsigned l,unsigned r,unsigned d){
	if( r<=a || b<=l ) return;
	if( a<=l && r<=b ){
		if( type == 'O' ) st[d].obstacle_tag += w;
		else st[d].shape_tag += w;
		return;
	}
	unsigned mid = (l+r)/2;
	insert(a,b,w,type,l,mid,d*2);
	insert(a,b,w,type,mid,r,d*2+1);
	up(l,r,d);
}
void swape_line::set_deep(unsigned a,unsigned b,int dep,unsigned l,unsigned r,unsigned d){
	if( r <= a || b <= l )return;
	if( a <= l && r <= b ){
		st[d].deep_tag = dep;
		return;
	}
	unsigned mid = (l+r)/2;
	down(d);
	set_deep(a,b,dep,l,mid,d*2);
	set_deep(a,b,dep,mid,r,d*2+1);
	deep_up(d);
}
tuple<int,char,unsigned> swape_line::findL(unsigned x,char type,unsigned l,unsigned r,unsigned d){
	unsigned shape_sum = get_shape_sum(l,r,d);
	unsigned obstacle_sum = get_obstacle_sum(l,r,d);
	int deep = get_deep(d);
	if( type == 0 ){
		if( shape_sum == r-l ) return make_tuple(deep,'S',l);
		if( obstacle_sum == r-l ) return make_tuple(deep,'O',l);
		if( shape_sum == 0 && obstacle_sum == 0 ) return make_tuple(deep,'W',l);
		unsigned mid = (l+r)/2;
		down(d);
		if( x <= mid ) return findL(x,0,l,mid,d*2);
		auto R = findL(x,0,mid,r,d*2+1);
		if( get<2>(R) != mid ) return R;
		auto L = findL(mid,get<1>(R),l,mid,d*2);
		if( get<1>(L) == 'X' ) return R;
		return L;
	}else{
		if( type == 'S' ){
			if( obstacle_sum == r-l ) return make_tuple(-1,'X',0);
			if( shape_sum == 0 && obstacle_sum == 0 ) return make_tuple(-1,'X',0);
			if( shape_sum == r-l ) return make_tuple(deep,'S',l);
		}else if( type == 'O' ){
			if( shape_sum == r-l ) return make_tuple(-1,'X',0);
			if( shape_sum == 0 && obstacle_sum == 0 ) return make_tuple(-1,'X',0);
			if( obstacle_sum == r-l ) return make_tuple(deep,'O',l);
		}else{
			if( obstacle_sum == r-l ) return make_tuple(-1,'X',0);
			if( shape_sum == r-l ) return make_tuple(-1,'X',0);
			if( shape_sum == 0 && obstacle_sum == 0 ) return make_tuple(deep,'W',l);
		}
		down(d);
		unsigned mid = (l+r)/2;
		if( x <= mid ) return findL(x,type,l,mid,d*2);
		auto R = findL(x,type,mid,r,d*2+1);
		if( get<1>(R) == 'X' ) return R;
		if( get<2>(R) != mid ) return R;
		auto L = findL(mid,type,l,mid,d*2);
		if( get<1>(L) == 'X' ) return R;
		return L;
	}
}
tuple<int,char,unsigned> swape_line::findR(unsigned x,char type,unsigned l,unsigned r,unsigned d){
	unsigned shape_sum = get_shape_sum(l,r,d);
	unsigned obstacle_sum = get_obstacle_sum(l,r,d);
	int deep = get_deep(d);
	if( type == 0 ){
		if( shape_sum == r-l ) return make_tuple(deep,'S',r);
		if( obstacle_sum == r-l ) return make_tuple(deep,'O',r);
		if( shape_sum == 0 && obstacle_sum == 0 ) return make_tuple(deep,'W',r);
		unsigned mid = (l+r)/2;
		down(d);
		if( x >= mid ) return findR(x,0,mid,r,d*2+1);
		auto L = findR(x,0,l,mid,d*2);
		if( get<2>(L) != mid ) return L;
		auto R = findR(mid,get<1>(L),mid,r,d*2+1);
		if( get<1>(R) == 'X' ) return L;
		return R;
	}else{
		if( type == 'S' ){
			if( obstacle_sum == r-l ) return make_tuple(-1,'X',0);
			if( shape_sum == 0 && obstacle_sum == 0 ) return make_tuple(-1,'X',0);
			if( shape_sum == r-l ) return make_tuple(deep,'S',r);
		}else if( type == 'O' ){
			if( shape_sum == r-l ) return make_tuple(-1,'X',0);
			if( shape_sum == 0 && obstacle_sum == 0 ) return make_tuple(-1,'X',0);
			if( obstacle_sum == r-l ) return make_tuple(deep,'O',r);
		}else{
			if( obstacle_sum == r-l ) return make_tuple(-1,'X',0);
			if( shape_sum == r-l ) return make_tuple(-1,'X',0);
			if( shape_sum == 0 && obstacle_sum == 0 ) return make_tuple(deep,'W',r);
		}
		unsigned mid = (l+r)/2;
		down(d);
		if( x >= mid ) return findR(x,type,mid,r,d*2+1);
		auto L = findR(x,type,l,mid,d*2);
		if( get<1>(L) == 'X' ) return L;
		if( get<2>(L) != mid ) return L;
		auto R = findR(mid,type,mid,r,d*2+1);
		if( get<1>(R) == 'X' ) return L;
		return R;
	}
}
bool swape_line::is_same(unsigned a,unsigned b,char type,unsigned l,unsigned r,unsigned d){
	if( r<=a || b<=l ) return 1;
	unsigned shape_sum = get_shape_sum(l,r,d);
	unsigned obstacle_sum = get_obstacle_sum(l,r,d);
	if( shape_sum == r-l ) return type=='S';
	if( obstacle_sum == r-l ) return type=='O';
	if( shape_sum==0 && obstacle_sum==0 ) return 0;
	unsigned mid = (l+r)/2;
	return is_same(a,b,type,l,mid,d*2)&&is_same(a,b,type,mid,r,d*2+1);
}
void swape_line::find_seg(unsigned a,unsigned b,unsigned l,unsigned r,unsigned d){
	if( r<=a || b<=l ) return;
	unsigned shape_sum = get_shape_sum(l,r,d);
	unsigned obstacle_sum = get_obstacle_sum(l,r,d);
	int deep = get_deep(d);
	unsigned L = max(l,a), R = min(r,b); 
	if( shape_sum == r-l ){
		segments.push_back(make_tuple(L,R,'S',deep));
		return;
	}
	if( obstacle_sum == r-l ){
		segments.push_back(make_tuple(L,R,'O',deep));
		return;
	}
	if( shape_sum==0 && obstacle_sum==0 ){
		segments.push_back(make_tuple(L,R,'W',deep));
		return;
	}
	unsigned mid = (l+r)/2;
	down(d);
	find_seg(a,b,l,mid,d*2);
	find_seg(a,b,mid,r,d*2+1);
}
