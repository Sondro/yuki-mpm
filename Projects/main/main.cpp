#include "globals.h"
#include "gridData.h"
#include "sampler.h"
int main(int argc, char **argv) {
	T n = 3;
	T r = 0.5;
	T k = 3;
	Sampler<T> *s = new Sampler<T>(n, r, k);
	return 0;
}