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
	Sampler(int numTiles, T n, T r, T k) : numTiles(numTiles), n(n), r(r), k(k),
                                           inverseCellWidth(std::sqrt(n) / r),
                                           gridResolution(std::floor(std::sqrt(n) / r)) {
        int seed = int(n * r * k);
        std::srand(seed);
        generateMasterTile();
        for (int i = 0; i < numTiles; i++) {
            tileSet.push_back(generatePoissonDistr(masterTile));
        }
 	}
    bool validPointSet(std::vector<vec3> allPoints) {
        for (vec3 &p : allPoints) {
            for (vec3 &o : allPoints) {
                if (p != o) {
                    if (p[0] = -1) continue;
                    if (o[0] = -1) continue;
                    if ((p - o).norm() < r) {
                        std::cout << "Points too close!" << std::endl;
                        std::cout << (p - o).norm() << std::endl;
                        std::cout << r << std::endl;
                        return false;
                    }
                }
            }
        }
        return true;
    }

    bool testTileCubeDistribution() {
        // tile each of the grids
        Partio::ParticlesDataMutable *parts = Partio::create();
        Partio::ParticleAttribute posH;

        posH = parts->addAttribute("position", Partio::VECTOR, 3);
        std::string filename = "multipleTiles.bgeo";
        std::vector<vec3> allPoints;
        for (int n = 0; n < numTiles; n++) {
            std::vector<vec3> data = tileSet[n];
            for (int i = 0; i < data.size(); ++i) {
                if (data[i][0] == -1) {
                    continue;
                }

                int idx = parts->addParticle();
                float *p = parts->dataWrite<float>(posH, idx);
                for (int k = 0; k < 3; ++k) {
                    p[k] = data[i][k];
                    switch (n) {
                        case 0:
                            break;
                        case 1:
                            if (k == 0) {
                                p[k] = 1 + data[i][k];
                            }
                            break;
                        case 2:
                            if (k == 1) {
                                p[k] = 1 + data[i][k];
                            }
                            break;
                        case 3:
                            if (k == 0 || k == 1) {
                                p[k] = 1 + data[i][k];
                            }
                            break;
                    }
                }
                allPoints.push_back(vec3(p[0], p[1], p[2]));
            }

        }

        Partio::write(filename.c_str(), *parts);
        parts->release();

        return validPointSet(allPoints);
    }
    bool testTileRowDistribution() {
        // tile each of the grids
        Partio::ParticlesDataMutable *parts = Partio::create();
        Partio::ParticleAttribute posH;

        posH = parts->addAttribute("position", Partio::VECTOR, 3);
        std::string filename = "multipleTiles.bgeo";
        std::vector<vec3> allPoints;
        for (int n = 0; n < numTiles; n++) {
            std::vector<vec3> data = tileSet[n];
            for (int i = 0; i < data.size(); ++i) {
                if (data[i][0] == -1) {
                    continue;
                }

                int idx = parts->addParticle();
                float *p = parts->dataWrite<float>(posH, idx);
                for (int k = 0; k < 3; ++k) {
                    p[k] = k == 0? data[i][k] + n : data[i][k];
                }
                allPoints.push_back(vec3(p[0], p[1], p[2]));
            }

        }

        Partio::write(filename.c_str(), *parts);
        parts->release();

        return validPointSet(allPoints);
    }
 	std::vector<vec3> generatePoissonDistr() {
        // initialize 3D grid
        std::vector<vec3> gridData(gridResolution * gridResolution * gridResolution, vec3(-1, -1, -1));
        std::vector<vec3> activeSamples;

        // initialize rng sampler
        vec3 x0(rng(), rng(), rng());

        std::vector<vec3> initialSamples;
        initialSamples.push_back(x0);
        return generatePoissonDistr(initialSamples);
 	}

    std::vector<vec3> generatePoissonDistr(std::vector<vec3> initialSamples) {
        // initialize 3D grid
        std::vector<vec3> gridData(gridResolution * gridResolution * gridResolution, vec3(-1, -1, -1));

        // add initial samples to grid
        for (vec3 sample : initialSamples) {
            int idx = getSampleIndex(sample);
            gridData[idx] = sample;
        }

        // add initial samples to active samples
        std::vector<vec3> activeSamples(initialSamples);

        while (activeSamples.size() > 0) {
            int idx = rng() * activeSamples.size();
            vec3 xi = activeSamples[idx];
            bool found = false;
            for (int i = 0; i < k; ++i) {
                vec3 sample = generateNeighborSample(xi, r);
                if (inGrid(sample) && isFarEnough(gridData, sample, r)) {
                    int newIdx = getSampleIndex(sample);
                    if (newIdx > 0 && newIdx < gridData.size()) {
                        found = true;
                        gridData[newIdx] = sample;
                        //printf("%f, %f, %f\n", sample[0], sample[1], sample[2]);
                        numSamples++;
                        activeSamples.push_back(sample);
                    }
                }
            }
            if (!found) {
                activeSamples.erase(activeSamples.begin() + idx);
            }
        }
        return gridData;
    }

    void generateMasterTile() {
        std::vector<vec3> gridData = generatePoissonDistr();
        // create master tile by getting list of points closer to edge than disk radius
        // test x boundaries
        for (int y = 0; y < gridResolution; y++) {
            for (int z = 0; z < gridResolution; z++) {
                // test lower x bound
                int idx = gridIndex3Dto1D(0, y, z);
                vec3 pt = gridData[idx];
                vec3 closestBoundaryPt(0, pt[1], pt[2]);
                T distanceToBoundary = (closestBoundaryPt - pt).norm();
                if (distanceToBoundary < r) {
                    masterTile.push_back(pt);
                }

                // test upper x bound
                idx = gridIndex3Dto1D(gridResolution - 1, y, z);
                pt = gridData[idx];
                closestBoundaryPt = vec3(1, pt[1], pt[2]);
                distanceToBoundary = (closestBoundaryPt - pt).norm();
                if (distanceToBoundary < r) {
                    masterTile.push_back(pt);
                }
            }
        }


        // test y boundaries
        for (int x = 0; x < gridResolution; x++) {
            for (int z = 0; z < gridResolution; z++) {
                int idx = gridIndex3Dto1D(x, 0, z);
                vec3 pt = gridData[idx];
                vec3 closestBoundaryPt(pt[0], 0, pt[2]);
                T distanceToBoundary = (closestBoundaryPt - pt).norm();
                if (distanceToBoundary < r) {
                    masterTile.push_back(pt);
                }

                idx = gridIndex3Dto1D(x, gridResolution - 1, z);
                pt = gridData[idx];
                closestBoundaryPt = vec3(pt[0], 1, pt[2]);
                distanceToBoundary = (closestBoundaryPt - pt).norm();
                if (distanceToBoundary < r) {
                    masterTile.push_back(pt);
                }
            }
        }

        // test z boundaries
        for (int x = 0; x < gridResolution; x++) {
            for (int y = 0; y < gridResolution; y++) {
                int idx = gridIndex3Dto1D(x, y, 0);
                vec3 pt = gridData[idx];
                vec3 closestBoundaryPt(pt[0], pt[1], 0);
                T distanceToBoundary = (closestBoundaryPt - pt).norm();
                if (distanceToBoundary < r) {
                    masterTile.push_back(pt);
                }

                idx = gridIndex3Dto1D(x, y, gridResolution - 1);
                pt = gridData[idx];
                closestBoundaryPt = vec3(pt[0], pt[1], 1);
                distanceToBoundary = (closestBoundaryPt - pt).norm();
                if (distanceToBoundary < r) {
                    masterTile.push_back(pt);
                }
            }
        }
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
 	bool isFarEnough(std::vector<vec3> data, vec3 sample, T r) {
        int sampleIdx = getSampleIndex(sample);
 		vec3 gridPos = sample * inverseCellWidth;
 		vec3 gridIdx = floor(gridPos);

        int search = 2;
 		for (int z = -search; z <= search; z++) {
			for (int y = -search; y <= search; y++) {
				for (int x = -search; x <= search; x++) {
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

 	vec3 getTiled(vec3 p) {
 		vec3 v = p;
 		for (int i = 0; i < 3; i++) {
 			if (v[i] < 0) v[i] += 1;
 			else if(v[i] > 1) v[i] -= 1;
 		}
 		
 		return v;
 	}
    void saveMasterTile(std::string filename) {
        Partio::ParticlesDataMutable *parts = Partio::create();
        Partio::ParticleAttribute posH;

        posH = parts->addAttribute("position", Partio::VECTOR, 3);

        for (int i = 0; i < masterTile.size(); ++i) {
            int idx = parts->addParticle();
            float *p = parts->dataWrite<float>(posH, idx);
            for (int k = 0; k < 3; ++k) {
                p[k] = masterTile[i][k];
            }
        }

        Partio::write(filename.c_str(), *parts);
        parts->release();
    }
 	// write particles to BGEO files
 	void saveSamples(std::vector<vec3> data, std::string filename) {
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

    std::vector<vec3> masterTile;
    std::vector<std::vector<vec3>> tileSet;
	T gridResolution;
	T inverseCellWidth;
    T n, r, k;
    int numSamples = 0;
    int numTiles;
};