#include "DataLoader.h"

int main()
{
    DataSet d;
    std::ifstream fin("a.in");
    d.load( fin );
}
