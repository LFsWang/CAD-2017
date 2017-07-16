struct statementP1{
	u32 a,b1,b2;//b2 > b1
	s32 inout;//-1 or 1
	statementP1(u32 a,u32 b1,u32 b2,s32 inout):a(a),b1(b1),b2(b2),inout(inout){}
	bool operator<(const statementP1 &B)const
	{
		return std::make_tuple(a,inout) < std::make_tuple(B.a,B.inout);
	}
};

struct Statemant_2D_VG{
	s32 type;
	u32 a; // x
	u32 b1;// y
	u32 b2;
	//1: out line
	//2: point
	//3: in line
	Statemant_2D_VG(s32 type,u32 a,u32 b1,u32 b2=0):type(type),a(a),b1(b1),b2(b2){}
	bool operator<(const Statemant_2D_VG &B)const{
		return std::make_tuple(a,type)<std::make_tuple(B.a,B.type);
	}
};