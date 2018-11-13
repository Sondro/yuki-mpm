#include "globals.h"
#include "gridData.h"
#include "sampler.h"
#include "samplertest.h"
#include "macGrid.h"
#include "simulation.h"

int main(int argc, char **argv) {
    T n = 3;
	T r2 = 1.0 / 2.0;
	T k = 30;

    // set up psuedo random sampling
    int seed = int(n * r2 * k);
    std::srand(seed);

	Sampler<T> *s = new Sampler<T>(n, r2 / 2.0, k);

    vec3i dim(10, 10, 10);
    vec3 initialPos;
    initialPos << dim[0] / 2, 7, dim[2] / 2;
    SamplerGrid<T> result = s->generatePoissonDistr();
    PointList samples = result.getAllSamples();
    s->saveSamples(samples, "samples.bgeo");
    std::vector<Particle<T>> particles;
    for (const auto &sample : samples) {
        particles.push_back(Particle<T>(sample));
    }
    Simulation sim(dim[0], dim[1], dim[2], initialPos, particles);
    std::cout << "running" << std::endl;
    sim.run();
	return 0;
}
