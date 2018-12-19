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
    int test = 0;
    // set up psuedo random sampling
    int seed = int(n * r2 * k);
    std::srand(seed);

    bool loadSamples = true;
	Sampler<T> *s;
	PointList samples;
	if (loadSamples) {
		// Loading samples
		s = new Sampler<T>(n, r2 / 2.0, k, "poisson_unit_cube_9443samples.bgeo");
		samples = s->unitCube.getAllSamples();
		for (auto &sample : samples) {
			sample /= 2;
		}
	}
	if (!loadSamples) {
		// Generating samples
		s = new Sampler<T>(n, r2 / 2.0, k);
		samples = s->unitCube.getAllSamples();
	}

    vec3i dim(X_CELL_COUNT, Y_CELL_COUNT, Z_CELL_COUNT);
    vec3 initialPos;
    initialPos << 1.5, 5.5, 1.5;
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
    sim.run();
	return 0;
}
