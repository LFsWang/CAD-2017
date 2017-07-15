#pragma once

struct statementP1{
	u32 a,b1,b2;//b2 > b1
	s32 inout;//-1 or 1
	statementP1(u32 a,u32 b1,u32 b2,s32 inout):a(a),b1(b1),b2(b2),inout(inout){}
	bool operator<(const statementP1 &B)const
	{
		return std::make_tuple(a,inout) < std::make_tuple(B.a,B.inout);
	}
};

struct swape_line_P1{
	struct node{
		u32 tag,area;
		node():tag(0),area(0){};
	};
	
	u32 n;
	std::vector<node> st;
	std::vector<std::pair<u32,u32>> segments;
	//(l,r)
	
	u32 get_area_sum(u32 l,u32 r,u32 d);
	void up(u32 l,u32 r,u32 d);
	
	void init(u32 _n);
	
	void insert(u32 a,u32 b,s32 w,u32 l,u32 r,u32 d);
	void find_seg(u32 a,u32 b,u32 l,u32 r,u32 d);
	
	void set_seg();
};

u32 swape_line_P1::get_area_sum(u32 l,u32 r,u32 d)
{
	if(st[d].tag) return r - l;
	return st[d].area;
}
void swape_line_P1::up(u32 l,u32 r,u32 d)
{
	u32 mid = (l+r)/2;
	st[d].area = get_area_sum(l,mid,d*2) + get_area_sum(mid,r,d*2+1);
}
void swape_line_P1::init(u32 _n)
{
	st=std::vector<node>( (n = _n) * 4 );
	segments.clear();
}
void swape_line_P1::insert(u32 a,u32 b,s32 w,u32 l,u32 r,u32 d)
{
	if( r<=a || b<=l ) return;
	if( a<=l && r<=b )
	{
		st[d].tag += w;
		return;
	}
	u32 mid = (l+r)/2;
	insert(a,b,w,l,mid,d*2);
	insert(a,b,w,mid,r,d*2+1);
	up(l,r,d);
}
void swape_line_P1::find_seg(u32 a,u32 b,u32 l,u32 r,u32 d)
{
	if( r<=a || b<=l ) return;
	u32 area = get_area_sum(l,r,d);
	u32 L = std::max(l,a), R = std::min(r,b); 
	if( area == r-l ) return;
	if( area==0 )
	{
		segments.emplace_back(L,R);
		return;
	}
	u32 mid = (l+r)/2;
	find_seg(a,b,l,mid,d*2);
	find_seg(a,b,mid,r,d*2+1);
}

void swape_line_P1::set_seg()
{
	std::vector<std::pair<u32,u32>> tmp;
	for(auto &i:segments)
	{
		if(tmp.size()&&tmp.back().second==i.first)
		{
			tmp.back().second=i.second;
		}
		else
		{
			tmp.emplace_back(i);
		}
	}
	segments.swap(tmp);
}