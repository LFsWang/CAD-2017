#include<iostream>
#include<fstream>
#include<cstdint>
#include<sstream>
#include<vector>

#include "DataLoader.h"

template<typename Ta,typename Tb>
static void range_to_value(std::string str,Ta &a,Tb &b)
{
    str = str.substr(1);
    str.pop_back();
    for(char &c:str)
        if(c==',')c=' ';
    std::stringstream ss(str);
    ss>>a>>b;
}

std::istream& operator>>(std::istream& in,DataSet::point &p)
{
    std::string buf;
    in>>buf;
    range_to_value(buf,p.x,p.y);
    return in;
}

std::ostream& operator<<(std::ostream& out,DataSet::point &p)
{
    out<<'('<<p.x<<','<<p.y<<')';
    return out;
}

static int rm_first_char(std::string s)
{
    return std::stoi(s.substr(1));
}

void DataSet::load(std::istream &fin)
{
    std::string buf;
    std::stringstream ss;

    auto getline_to_ss = [&](){
        std::getline(fin,buf);
        ss.clear();
        ss.str(buf);
    };

    getline_to_ss(); //ViaCost
    ss>>buf>>buf>>viacost;
    
    getline_to_ss(); //Spacing
    ss>>buf>>buf>>spacing;
    
    getline_to_ss(); //Boundary
    ss>>buf>>buf>>boundary.first>>boundary.second; 

    getline_to_ss(); //MetalLayers
    ss>>buf>>buf>>metal_layers;
    for(s32 i=1;i<=metal_layers;++i)
    {
        Obstacles[i].clear();
        RoutedShape[i].clear();
        RoutedVia[i].clear();
    }
    
    getline_to_ss();
    ss>>buf>>buf>>routed_shapes;
    
    getline_to_ss(); //RoutedVias
    ss>>buf>>buf>>routed_vias;

    getline_to_ss(); //Obstacles  
    ss>>buf>>buf>>obstacles;

    DataSet::point tmp_point;
    int lay;
    for(s32 i=0;i<routed_shapes;++i)
    {
        getline_to_ss();
        ss>>buf>>buf;
        lay = rm_first_char(buf);
        RoutedShape[lay].emplace_back(tmp_point,tmp_point);
        ss>>RoutedShape[lay].back().first>>RoutedShape[lay].back().second;
    }
    
    for(s32 i=0;i<routed_vias;++i)
    {
        getline_to_ss();
        ss>>buf>>buf;
        lay = rm_first_char(buf);
        RoutedVia[lay].emplace_back(tmp_point);
        ss>>RoutedVia[lay].back();  
    }
    
    for(s32 i=0;i<obstacles;++i)
    {
        getline_to_ss();
        ss>>buf>>buf;
        lay = rm_first_char(buf);
        Obstacles[lay].emplace_back(tmp_point,tmp_point);
        ss>>Obstacles[lay].back().first>>Obstacles[lay].back().second;
    }
}


