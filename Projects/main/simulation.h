#pragma once

#include "globals.h"
#include "macGrid.h"

class Simulation {
public:
	Simulation(int i, int j, int k,
			   mat4 transform,
			   std::vector<Particle<T>> particle,
			   double seconds = 20.0,
			   double fps = 24.0);
	void run();

	MACGrid<T> grid;
	double seconds;
	double fps;
	double frames;

private:
};