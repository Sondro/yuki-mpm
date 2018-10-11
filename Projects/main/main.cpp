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
