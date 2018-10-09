#include "globals.h"
#include "gridData.h"
#include "sampler.h"
int main(int argc, char **argv) {
	T n = 3;
	T r = 0.15;
	T k = 30;
    int numTiles = 4;
	Sampler<T> *s = new Sampler<T>(numTiles, n, r, k);
    //s->saveMasterTile("poisson_masterTile.bgeo");
    //std::cout << s->validPointSet(s->generatePoissonDistr());
    std::cout << s->testTileCubeDistribution();
	return 0;
}