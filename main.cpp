#include "DataLoader.h"
#include "BuildVisingGraph.h"

int main()
{
    DataSet d;
    std::ifstream fin("a.in");
    d.load( fin );
	d.set_spacing_on_Obstacles();
	VisingGraph v;
	v.build(d);
}
