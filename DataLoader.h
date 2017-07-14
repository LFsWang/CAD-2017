#pragma once
#include<iostream>
#include<fstream>
#include<cstdint>
#include<sstream>
#include<vector>

using u64 = std::uint_fast64_t;
using s64 = std::int_fast64_t;
using u32 = std::uint_fast32_t;
using s32 = std::int_fast32_t;
using u16 = std::uint_fast16_t;
using s16 = std::int_fast16_t;
using u8 = std::uint_fast8_t;
using s8 = std::int_fast8_t;
using std::cout;
using std::endl;

class DataSet{
public:
    static const int LIMIT_LAYER = 10+1;//for 1-base
    struct point{
        s64 x;
        s64 y;
    };
    s64 viacost;
    s64 spacing;
    std::pair<point,point> boundary;

    s32 metal_layers;
    s32 routed_shapes;
    s32 routed_vias;
    s32 obstacles;

    std::vector<std::pair<point,point>> Obstacles  [LIMIT_LAYER];
    std::vector<std::pair<point,point>> RoutedShape[LIMIT_LAYER];
    std::vector<point> RoutedVia  [LIMIT_LAYER];
       
    
    void load(std::istream &fin);
	void set_spacing_on_Obstacles();
};
std::istream& operator>>(std::istream& in,DataSet::point &p);
std::ostream& operator<<(std::ostream& out,DataSet::point &p);


