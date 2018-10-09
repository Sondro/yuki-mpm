#include "globals.h"
#include "gridData.h"
#include "sampler.h"
int main(int argc, char **argv) {
	T n = 3;
	T r = 0.1;
	T k = 30;
	Sampler<T> *s = new Sampler<T>(n, r, k);
	return 0;
}