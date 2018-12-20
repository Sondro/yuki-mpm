#include <iostream>
#include <iomanip>
#include <chrono>
#include "simulation.h"

Simulation::Simulation(int i, int j, int k, mat4 transform, std::vector<Particle<T>> particles, double seconds, double fps) :
	grid(i, j, k),
	seconds(seconds),
	fps(fps),
	frames(fps * seconds) {
    grid.particles = particles;
    for (auto &p : grid.particles) {
        vec4 tmp;
        tmp << p.pos[0], p.pos[1], p.pos[2], 1.0;
        tmp = transform * tmp;
        p.pos[0] = tmp[0];
        p.pos[1] = tmp[1];
        p.pos[2] = tmp[2];
    }
}

namespace {
void printTime(float seconds) {
	int i = static_cast<int>(std::round(seconds));
	int h = i / 3600;
	int m = (i % 3600) / 60;
	int s = i % 60;
	if (h > 0) {
		std::cout << h << " h ";
	}
	if (h > 0 || m > 0) {
		std::cout << m << " m ";
	}
	std::cout << s << " s";
}
}

void Simulation::run() {
	std::cout << "Beginning simulation..." << std::endl;
	double cumTime = 0.0;
	double simTime = 0.0;
	for (int frame = 0; frame < frames; ++frame) {
		std::cout << "  Writing frame " << std::setw(3) << std::setfill('0') << frame + 1 <<
				  " at timestep " << cumTime << std::endl;
        grid.writeFrame();
        if (frame != 0) {
			double avg = simTime * 1e-6 / frame;
			std::cout << "    Estimated time remaining: ";
			printTime(avg * (frames - frame));
			std::cout << "\n";
		}
		auto start = std::chrono::high_resolution_clock::now();

		for (double timestep = 0; timestep < 1.0 / fps; timestep += dt, cumTime += dt) {
			grid.stepSimulation();
		}

		auto finish = std::chrono::high_resolution_clock::now();
		auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
		simTime += microseconds.count();
	}
}