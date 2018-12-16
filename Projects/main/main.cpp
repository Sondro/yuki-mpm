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

    vec3i dim(X_CELL_COUNT, Y_CELL_COUNT, Z_CELL_COUNT);
    vec3 initialPos;
    initialPos << 3, 5.5, 3;
    SamplerGrid<T> result = s->generatePoissonDistr();
    PointList samples = result.getAllSamples();
    s->saveSamples(samples, "samples.bgeo");
    std::vector<Particle<T>> particles;
    for (const auto &sample : samples) {
    	Particle<T> p(sample);
		p.vol = VOLUME / samples.size();
		p.mass = DENSITY * p.vol;
        particles.push_back(p);
    }
    Simulation sim(dim[0], dim[1], dim[2], initialPos, particles);
    sim.grid.drawGrid();
    std::cout << "running" << std::endl;
    sim.run();
	return 0;
}
