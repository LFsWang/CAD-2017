#pragma once
// 7122
struct BIT{
	s32 n;
	std::vector<s32> st;
	void init(s32 _n){
		st.resize((n=_n)+1);
	}
	void clear(){
		for(auto &i:st)
			i=0;
	}
	void add(s32 id,s32 data){
		for(;id<=n;id+=(id&-id))
			st[id]+=data;
	}
	s32 get_sum(s32 id){
		s32 res=0;
		for(;id;id-=(id&-id))
			res+=st[id];
		return res;
	}
};
