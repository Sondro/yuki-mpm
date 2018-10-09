template <typename T>
// Fast Poisson Disk Sampling in Arbitrary Dimensions, Robert Bridson
class Sampler {
public:
	/*
	 * param: T n
	 	dimension of the sampling domain
	 * param: T r
	 	minimum distance r between samples
	 * param: T k 
	 	limit of samples to choose before rejection 
	 	in algorithm
	 */
	Sampler(T n, T r, T k) { 
		inverseCellWidth = std::sqrt(n) / r;
		gridResolution = std::floor(std::sqrt(n) / r);
		// initialize 3D grid
		data = std::vector<vec3>(gridResolution * gridResolution * gridResolution, vec3(-1, -1, -1));
		
		// initialize rng sampler
		int seed = int(n * r * k);
		std::srand(seed);
		vec3 x0(rng(), rng(), rng());

		// add first sample to grid data
        int flatIdx = getSampleIndex(x0);
        std::cout << flatIdx << std::endl;
		data[flatIdx] = x0;

		// initialize active list
		activeSamples.push_back(x0);

		while (activeSamples.size() > 0) {
			int idx = rng() * activeSamples.size();
			vec3 xi = activeSamples[idx];
			bool found = false;
			for (int i = 0; i < k; ++i) {
				vec3 sample = generateNeighborSample(xi, r);
				if (inGrid(sample) && isFarEnough(sample, r)) {
					int newIdx = getSampleIndex(sample);
					if (newIdx > 0 && newIdx < data.size()) {
                        found = true;
						data[newIdx] = sample;
                        printf("%f, %f, %f\n", sample[0], sample[1], sample[2]);
                        numSamples++;
						activeSamples.push_back(sample);
					}
				}
			}
			if (!found) {
                //std::cout << "removed element" << std::endl;
                //std::cout << activeSamples.size() << std::endl;
				activeSamples.erase(activeSamples.begin() + idx);
			}
		}

        std::cout << "finished generating samples" << std::endl;
        std::cout <<  numSamples << std::endl;
 	}

 	vec3 generateNeighborSample(vec3 xi, T r) {
 		// randomly generate spherical coordinates
 		T radius = (rng() + 1) * r;
 		T angle1 = 2.0 * PI * rng();
 		T angle2 = 2.0 * PI * rng();
 		T x = xi[0] + radius * std::cos(angle1) * std::sin(angle2);
 		T y = xi[1] + radius * std::sin(angle1) * std::sin(angle2);
 		T z = xi[2] + radius * std::cos(angle2);
 		return vec3(x, y, z);
 	}

    int getSampleIndex(vec3 sample) {
        vec3 gridPos = sample * inverseCellWidth;
        vec3 gridIdx = floor(gridPos);

        int flatIdx = gridIndex3Dto1D(gridIdx[0], gridIdx[1], gridIdx[2]);
        return flatIdx;
    }
    bool inGrid(vec3 sample) {
        for (int i = 0; i < 3; i++) {
            if (sample[i] < 0 || sample[i] > 1) return false;
        }
        return true;
    }
 	bool isFarEnough(vec3 sample, T r) {
        int sampleIdx = getSampleIndex(sample);
 		vec3 gridPos = sample * inverseCellWidth;
 		vec3 gridIdx = floor(gridPos);

 		for (int z = -1; z <= 1; z++) {
			for (int y = -1; y <= 1; y++) {
				for (int x = -1; x <= 1; x++) {
					vec3 offset(x, y, z);
					vec3 neighborIdx = offset + gridIdx;
					if (neighborIdx[0] < 0 || neighborIdx[0] >= gridResolution ||
						neighborIdx[1] < 0 || neighborIdx[1] >= gridResolution ||
						neighborIdx[2] < 0 || neighborIdx[2] >= gridResolution) 
						continue;

					int flatNeighborIdx = gridIndex3Dto1D(neighborIdx[0], neighborIdx[1], neighborIdx[2]);
					if (data[flatNeighborIdx] == vec3(-1, -1, -1)) {

                    }
					else {
						vec3 neighborSample = data[flatNeighborIdx];
						vec3 v = neighborSample - sample;
						T distance = v.norm();
						if (distance < r) return false;
					}

				}
			}
		}
		return true;
 	}
 	int gridIndex3Dto1D(vec3 gridPos) {
 		int x = int(gridPos[0]);
 		int y = int(gridPos[1]);
 		int z = int(gridPos[2]);
 		return x + y * gridResolution + z * gridResolution * gridResolution;
 	}

 	int gridIndex3Dto1D(int x, int y, int z) {
 		return x + y * gridResolution + z * gridResolution * gridResolution;
 	}

 	// write particles to BGEO files
 	void saveSamples(std::string filename) {
 		Partio::ParticlesDataMutable *parts = Partio::create();
 		Partio::ParticleAttribute posH;

 		posH = parts->addAttribute("position", Partio::VECTOR, 3);

 		for (int i = 0; i < data.size(); ++i) {
 			if (data[i][0] == -1) {
 				continue;
 			}
 			int idx = parts->addParticle();
 			float *p = parts->dataWrite<float>(posH, idx);
 			for (int k = 0; k < 3; ++k) {
 				p[k] = data[i][k];
 			}
 		}

 		Partio::write(filename.c_str(), *parts);
 		parts->release();
 	}

 	T rng() { return ((T) std::rand() / (RAND_MAX)); }
	std::vector<vec3> data;
	std::vector<vec3> activeSamples;
	T gridResolution;
	T inverseCellWidth;
    int numSamples = 0;
};