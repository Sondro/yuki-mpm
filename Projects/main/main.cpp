#include "globals.h"
#include "sampler.h"
#include "samplertest.h"

int main(int argc, char **argv) {
	T n = 3;
	T r2 = 1.0 / 5.0;
	T k = 30;

    // set up psuedo random sampling
    int seed = int(n * r2 * k);
    std::srand(seed);

	Sampler<T> *s = new Sampler<T>(n, r2 / 2.0, k);
    std::cout << "Num samples: " << s->numSamples << std::endl;

    SamplerTest<T> test(s);
    test.testUnitPoissonCube();
    test.testTileCubeDistribution();
    test.testTileRowDistribution();
	return 0;
}


//using namespace voro;
//
//// This function returns a random floating point number between 0 and 1
//double rnd() {return double(rand())/RAND_MAX;}
//
//int main() {
//	double x,y,z,rsq,r;
//	voronoicell v;
//
//	// Initialize the Voronoi cell to be a cube of side length 2, centered
//	// on the origin
//	v.init(-1,1,-1,1,-1,1);
//
//	// Cut the cell by 250 random planes which are all a distance 1 away
//	// from the origin, to make an approximation to a sphere
//	for(int i=0;i<250;i++) {
//		x=2*rnd()-1;
//		y=2*rnd()-1;
//		z=2*rnd()-1;
//		rsq=x*x+y*y+z*z;
//		if(rsq>0.01&&rsq<1) {
//			r=1/sqrt(rsq);x*=r;y*=r;z*=r;
//			v.plane(x,y,z,1);
//		}
//	}
//
//	// Output the Voronoi cell to a file, in the gnuplot format
//	v.draw_gnuplot(0,0,0,"single_cell.gnu");
//}

