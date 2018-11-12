#include "globals.h"
#include "gridData.h"
#include "sampler.h"
#include "samplertest.h"
#include "macGrid.h"
#include "simulation.h"

int main(int argc, char **argv) {
    Simulation sim(vec3i(10, 10, 10));
    sim.run();
	return 0;
}
