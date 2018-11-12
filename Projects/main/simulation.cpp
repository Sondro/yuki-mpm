#include <iostream>
#include <iomanip>
#include <chrono>
#include "simulation.h"

Simulation::Simulation(vec3i gridDims, double seconds, double fps) :
	grid(gridDims),
	seconds(seconds),
	fps(fps),
	frames(fps * seconds) {}

void Simulation::run() {
	std::cout << "Beginning simulation..." << std::endl;
	double cumTime = 0.0;
	double simTime = 0.0;
	for (int frame = 0; frame < frames; ++frame) {
		std::cout << "  Writing frame " << std::setw(3) << std::setfill('0') << frame + 1 <<
				  " at timestep " << cumTime << std::endl;
		if (frame != 0) {
			double avg = simTime * 1e-6 / frame;
			std::cout << "    Estimated time remaining: " << avg * (frames - frame) << " s" << std::endl;
		}
		auto start = std::chrono::high_resolution_clock::now();
		grid.writeFrame();

		for (double timestep = 0; timestep < 1.0 / fps; timestep += dt, cumTime += dt) {
			grid.stepSimulation();
		}

		auto finish = std::chrono::high_resolution_clock::now();
		auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
		simTime += microseconds.count();
	}
}