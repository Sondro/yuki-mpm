#include "globals.h"
#include "sampler.h"
#include "samplertest.h"
int main(int argc, char **argv) {
	T n = 3;
	T r2 = 1.0 / 5.0;
	T k = 30;
    int numTiles = 4;
	Sampler<T> *s = new Sampler<T>(numTiles, n, r2 / 2.0, k);
    std::cout << "Num samples: " << s->numSamples << std::endl;

    SamplerTest<T> test(s);
    test.testTileCubeDistribution();
	return 0;
}