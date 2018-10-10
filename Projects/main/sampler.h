#pragma once

#include "globals.h"
#define DEBUG 0
#define GETX 1
#define GETY 1
#define GETZ 1

template <typename T>
class SamplerGrid {
    using vector = std::vector<PointList>;
public:
    SamplerGrid(int gridResolution, T cellWidth) : gridResolution(gridResolution),
    cellWidth(cellWidth), inverseCellWidth(1.0 / cellWidth) {
        mData = vector(gridResolution * gridResolution * gridResolution, PointList());
    }

    ~SamplerGrid() {}


    int flat(int x, int y, int z) const {
        return x + y * gridResolution + z * gridResolution * gridResolution;
    }

    virtual PointList &operator()(int i, int j, int k) {
        int idx = 0;
        idx = flat(i, j, k);
        return mData[idx];
    }

    virtual PointList operator()(int i, int j, int k) const {
        int idx = 0;
        idx = flat(i, j, k);
        return mData[idx];
    }

    virtual PointList &operator[](int idx) {
        return mData[idx];
    }

    virtual PointList operator[](int idx) const {
        return mData[idx];
    }


    int getCellIndex(vec3 sample) {
        vec3 gridPos = sample * inverseCellWidth;
        return flat(int(gridPos[0]), int(gridPos[1]), int(gridPos[2]));
    }

    bool addSample(vec3 pt) {
        int idx = getCellIndex(pt);
        if (idx < 0 || idx >= mData.size()) return false;
        (*this)[idx].push_back(pt);
        return true;
    }

    PointList getAllSamples() {
        PointList pts;
        for (int i = 0; i < gridResolution; i++) {
            for (int j = 0; j < gridResolution; j++) {
                for (int k = 0; k < gridResolution; k++) {
                    PointList pointList = (*this)(i, j, k);
                    for (vec3 &p : pointList) {
                        pts.push_back(p);
                    }
                }
            }
        }
        return pts;
    }

    int gridResolution;
    T cellWidth;
    T inverseCellWidth;
    vector mData;
};

template <typename T>
class Sampler {
public:
	/*
	 * Creates the master tile by taking the outer samples of an initial
	 * poisson distribution. Then creates n poisson unit cube distributions
	 * with the master tile as the initial active sample set.
	 *
	 * param: int numTiles
	 *  number of stampable poisson cubes to generate
	 * param: T n
	 *  dimension of the sampling domain
	 * param: T r
	 *  minimum distance r between samples
	 * param: T k 
	 *  limit of samples to choose before rejection in algorithm
	 */
	Sampler(int numTiles, T n, T r, T k) : numTiles(numTiles), n(n), r(r), k(k),
                                           cellWidth(2.0 * r),
                                           inverseCellWidth(1.0 / (2.0 * r)),
                                           gridResolution(1.0 / (2.0 * r)) {
        // set up psuedo random sampling
        int seed = int(n * r * k);
        std::srand(seed);

//        SamplerGrid<T> result = generatePoissonDistr();
//        PointList samples = result.getAllSamples();
//        validPointSet(samples);
//        saveSamples(samples, "fixedPoisson.bgeo");

        // master tile is the outer samples for each cube
        generateMasterTile();
        saveSamples(masterTile, "masterTile.bgeo");

        // create poisson unit cubes
        for (int i = 0; i < numTiles; i++) {
            SamplerGrid<T> cube = generatePoissonDistr(masterTile);
            tileSet.push_back(cube.getAllSamples());
        }
 	}

 	SamplerGrid<T> generatePoissonDistr() {
        // initialize 3D grid
        std::vector<vec3> activeSamples;

        // initialize rng sampler
        vec3 x0(rng(), rng(), rng());

        std::vector<vec3> initialSamples;
        initialSamples.push_back(x0);
        return generatePoissonDistr(initialSamples);
 	}

    SamplerGrid<T> generatePoissonDistr(std::vector<vec3> initialSamples) {
        // initialize 3D grid
        SamplerGrid<T> gridData(gridResolution, cellWidth);

        // add initial samples to grid
        for (vec3 sample : initialSamples) {
            gridData.addSample(sample);
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
                    bool added = gridData.addSample(sample);
                    if (added) {
                        found = true;
                        activeSamples.push_back(sample);
                        numSamples++;
#if DEBUG
                        printf("%f, %f, %f\n", sample[0], sample[1], sample[2]);
#endif
                    }
                }
            }
            if (!found) {
                activeSamples.erase(activeSamples.begin() + idx);
#if DEBUG
                std::cout << "Removed sample. Active sample count: " << activeSamples.size() << std::endl;
#endif
            }
        }
        return gridData;
    }

    bool validPointSet(PointList allPoints) {
        bool testPassed = true;
        for (vec3 &p : allPoints) {
            for (vec3 &o : allPoints) {
                if (p != o) {
                    // return false if the points are closer than the poisson
                    // disk radius
                    if ((p - o).norm() < r) {
                        std::cout << "Points too close!" << std::endl;
                        std::cout << "Distance: " << (p - o).norm() << std::endl;
                        std::cout << "Poisson Radius: " << r << std::endl;

                        testPassed = false;
                    }
                }
            }
        }
        std::cout << "Points meet poisson disk radius constraint." << std::endl;
        return testPassed;
    }


    void generateMasterTile() {
        SamplerGrid<T> gridData = generatePoissonDistr();
        PointList ptList = gridData.getAllSamples();
#ifdef DEBUG
        std::cout << "Generating master tile..." << std::endl;
        std::string result = validPointSet(ptList) ? "Valid distr" : "Invalid distr.";
        std::cout << result << std::endl;
#endif
        // create master tile by getting list of points closer to edge than disk radius
        // test x boundaries
#if GETX
        for (int y = 0; y < gridResolution; y++) {
            for (int z = 0; z < gridResolution; z++) {
                // test lower x bounds
                PointList pts = gridData(0, y, z);
                for (vec3 &pt : pts) {
                    vec3 closestBoundaryPt(0, pt[1], pt[2]);
                    T distanceToBoundary = (closestBoundaryPt - pt).norm();
                    if (distanceToBoundary < r) {
                        masterTile.push_back(pt);
                    }
                }

                // test upper x bounds
                pts = gridData(gridResolution - 1, y, z);
                for (vec3 &pt : pts) {
                    vec3 closestBoundaryPt(1, pt[1], pt[2]);
                    T distanceToBoundary = (closestBoundaryPt - pt).norm();
                    if (distanceToBoundary < r) {
                        masterTile.push_back(pt);
                    }
                }

            }
        }
#endif

#if GETY
        // test y boundaries
        for (int x = 0; x < gridResolution; x++) {
            for (int z = 0; z < gridResolution; z++) {
                // test lower y bounds
                PointList pts = gridData(x, 0, z);
                for (vec3 &pt : pts) {
                    vec3 closestBoundaryPt(pt[0], 0, pt[2]);
                    T distanceToBoundary = (closestBoundaryPt - pt).norm();
                    if (distanceToBoundary < r) {
                        masterTile.push_back(pt);
                    }
                }

                // test upper y bounds
                pts = gridData(x, gridResolution - 1, z);
                for (vec3 &pt : pts) {
                    vec3 closestBoundaryPt(pt[0], 1, pt[2]);
                    T distanceToBoundary = (closestBoundaryPt - pt).norm();
                    if (distanceToBoundary < r) {
                        masterTile.push_back(pt);
                    }
                }
            }
        }
#endif

#if GETZ
        // test z boundaries
        for (int x = 0; x < gridResolution; x++) {
            for (int y = 0; y < gridResolution; y++) {
                // test lower y bounds
                PointList pts = gridData(x, y, 0);
                for (vec3 &pt : pts) {
                    vec3 closestBoundaryPt(pt[0], pt[1], 0);
                    T distanceToBoundary = (closestBoundaryPt - pt).norm();
                    if (distanceToBoundary < r) {
                        masterTile.push_back(pt);
                    }
                }

                // test upper y bounds
                pts = gridData(x, y, gridResolution - 1);
                for (vec3 &pt : pts) {
                    vec3 closestBoundaryPt(pt[0], pt[1], 1);
                    T distanceToBoundary = (closestBoundaryPt - pt).norm();
                    if (distanceToBoundary < r) {
                        masterTile.push_back(pt);
                    }
                }
            }
        }
#endif

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

    bool inGrid(vec3 sample) {
        for (int i = 0; i < 3; i++) {
            if (sample[i] < 0 || sample[i] > 1) return false;
        }
        return true;
    }

 	bool isFarEnough(SamplerGrid<T> data, vec3 sample, T r) {
 		vec3 gridPos = sample * inverseCellWidth;
 		vec3 gridIdx = floor(gridPos);

        int search = 1;
 		for (int z = -search; z <= search; z++) {
			for (int y = -search; y <= search; y++) {
				for (int x = -search; x <= search; x++) {
					vec3 offset = vec3(x, y, z);
					vec3 neighborIdx = offset + gridIdx;
                    int i = neighborIdx[0];
                    int j = neighborIdx[1];
                    int k = neighborIdx[2];
					if (i < 0 || i >= gridResolution ||
						j < 0 || j >= gridResolution ||
						k < 0 || k >= gridResolution)
						continue;

                    PointList cellPts = data(i, j, k);
                    if (!cellPts.empty()) {
                        for (vec3 &neighborSample : cellPts) {
                            vec3 v = neighborSample - sample;
                            T distance = v.norm();
                            if (distance < r) return false;
                        }
                    }
				}
			}
		}
		return true;
 	}

 	// write particles to BGEO files
 	void saveSamples(std::vector<vec3> data, std::string filename) {
 		Partio::ParticlesDataMutable *parts = Partio::create();
 		Partio::ParticleAttribute posH;

 		posH = parts->addAttribute("position", Partio::VECTOR, 3);

 		for (int i = 0; i < data.size(); ++i) {
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
	int gridResolution;
	T inverseCellWidth;
    T cellWidth;
    T n, r, k;
    int numSamples = 0;
    int numTiles;
};