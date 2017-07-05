#include"rebuild_input.h"
using namespace std;
input_data input;
int main(){
	input.init_by_stdin();
	auto output = rebuild_input_by_swape_line::rebuild_input(input);
	output.print();
	return 0;
}
/*
ViaCost = Cv
Spacing = S
Boundary = (LLx,LLy) (URx,URy)
#MetalLayers = W
#RoutedShapes = X
#RoutedVias = Y
#Obstacles = Z
RoutedShape Layer (LLx,LLy) (URx,URy)
...
RoutedVia Layer (x,y)
...
Obstacle Layer (LLx,LLy) (URx,URy)
*/

/*
ViaCost = 20
Spacing = 5
Boundary = (0,0) (1000,1000)
#MetalLayers = 2
#RoutedShapes = 7
#RoutedVias = 1
#Obstacles = 3
RoutedShape M1 (50,100) (250,150)
RoutedShape M1 (600,20) (750,140)
RoutedShape M1 (50,850) (250,900)
RoutedShape M1 (10,800) (500,995)
RoutedShape M2 (75,20) (200,750)
RoutedShape M2 (375,100) (575,600)
RoutedShape M2 (475,20) (670,450)
RoutedVia V1 (175,125)
Obstacle M1 (350,300) (650,750)
Obstacle M1 (50,350) (650,650)
Obstacle M2 (350,700) (950,800)
*/
