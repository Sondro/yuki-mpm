#pragma once

#include "globals.h"
#define DEBUG 1
#define GETX 1
#define GETY 1
#define GETZ 1

template <typename T>

/**
 * Uniform grid data structure to contain the samples.
 * @tparam T
 */
class SamplerGrid {
    using FlatGrid = std::vector<PointList>;
public:

    /**
     * Initialize the grid as a 2D vector for points.
     * @param gridResolution number of cells per dimension
     * @param cellWidth size of each cell
     */
    SamplerGrid(int gridResolution, T cellWidth) :
            gridResolution(gridResolution),
            cellWidth(cellWidth),
            inverseCellWidth(1.0 / cellWidth) {
        mData = FlatGrid(gridResolution * gridResolution * gridResolution,
                         PointList());
    }

    ~SamplerGrid() {}

    /**
     * Return the actual index in the mData array using the
     * 3D coordinates.
     */
    int flat(int x, int y, int z) const {
        return x + y * gridResolution + z * gridResolution * gridResolution;
    }

    /**
     * Get the points in the cell at i, j, k. i, j, k are indices.
     * @param i
     * @param j
     * @param k
     * @return list of samples contained in the cell at i, j, k.
     */
    virtual PointList &operator()(int i, int j, int k) {
        int idx = 0;
        idx = flat(i, j, k);
        return mData[idx];
    }

    /**
     * Get the points in the cell at i, j, k. i, j, k are indices.
     * @param i
     * @param j
     * @param k
     * @return list of samples contained in the cell at i, j, k.
     */
    virtual PointList operator()(int i, int j, int k) const {
        int idx = 0;
        idx = flat(i, j, k);
        return mData[idx];
    }

    /**
     * Get the list of points in the cell at idx.
     * @param idx
     * @return list of samples contained in mData at idx.
     */
    virtual PointList &operator[](int idx) {
        return mData[idx];
    }

    /**
     * Get the list of points in the cell at idx.
     * @param idx
     * @return list of samples contained in mData at idx.
     */
    virtual PointList operator[](int idx) const {
        return mData[idx];
    }

    /**
     * Get the flat index for the cell that this sample belongs to.
     * @param sample
     * @return flat index for the sample.
     */
    int getCellIndex(vec3 sample) {
        vec3 gridPos = sample * inverseCellWidth;
        return flat(int(gridPos[0]), int(gridPos[1]), int(gridPos[2]));
    }

    /**
     * Add the sample to the grid at the cell.
     * @param pt
     * @return whether the sample is successfully added.
     */
    bool addSample(vec3 pt) {
        int idx = getCellIndex(pt);
        if (idx < 0 || idx >= mData.size()) return false;
        (*this)[idx].push_back(pt);
        return true;
    }

    /**
     * Return a list of samples contained in this grid data structure.
     * @return lists of samples.
     */
    PointList getAllSamples() {
        PointList pts;
        for (int i = 0; i < gridResolution; i++) {
            for (int j = 0; j < gridResolution; j++) {
                for (int k = 0; k < gridResolution; k++) {
                    PointList pointList = (*this)(i, j, k);
                    for (vec3 p : pointList) {
                        pts.push_back(p);
                    }
                }
            }
        }
        return pts;
    }

    int gridResolution; // number of cells for each side.
    T cellWidth; // width of the cells in the grid.
    T inverseCellWidth; // 1 / cellWidth.

private:
    FlatGrid mData; // array of point lists. there is a list of points per cell.
};

/**
 * Class that produces a tile-able poisson distribution.
 * @tparam T
 */
template <typename T>
class Sampler {
public:
	/**
	 * Creates the master tile by taking the outer samples of an initial
	 * poisson distribution. Then creates n poisson unit cube distributions
	 * with the master tile as the initial active sample set.
	 */
	Sampler(int numTiles, T n, T r, T k) : numTiles(numTiles),
                                           n(n), r(r), k(k),
                                           cellWidth(2.0 * r),
                                           inverseCellWidth(1.0 / (2.0 * r)),
                                           gridResolution(1.0 / (2.0 * r)) {
        // set up psuedo random sampling
        int seed = int(n * r * k);
        std::srand(seed);

        // master tile is the outer samples for each cube
        generateMasterTile();
        saveSamples(masterTile, "masterTile.bgeo");

        // create poisson unit cubes
        for (int i = 0; i < numTiles; i++) {
            SamplerGrid<T> cube = generatePoissonDistr(masterTile);
            tileSet.push_back(cube.getAllSamples());
        }
 	}

    /**
     * Generate poisson samples using Fast Poisson Disk Sampling. Grid is
     * initialized with a single sample.
     * @return a sampler grid data structure filled with poisson disk samples.
     */
 	SamplerGrid<T> generatePoissonDistr() {
        // initialize 3D grid
        std::vector<vec3> activeSamples;

        // initialize rng sampler
        vec3 x0(rng(), rng(), rng());

        std::vector<vec3> initialSamples;
        initialSamples.push_back(x0);
        return generatePoissonDistr(initialSamples);
 	}

    /**
     * Generate poisson samples using Fast Poisson Disk Sampling. Grid is
     * initialized with a list of samples (master tile).
     * @param initialSamples list of samples already in grid before sampling begins
     * @return a sampler grid data structure filled with poisson disk samples.
     */
    SamplerGrid<T> generatePoissonDistr(std::vector<vec3> initialSamples) {
        // initialize 3D grid
        SamplerGrid<T> gridData(gridResolution, cellWidth);

        // add initial samples to grid
        for (vec3 sample : initialSamples) {
            gridData.addSample(sample);
        }

        // add initial samples to active samples
        std::vector<vec3> activeSamples(initialSamples);

        // continue to add samples to the grid around a sample until no more
        // samples that are far enough can be generated
        while (activeSamples.size() > 0) {
            // randomly choose a sample to generate other samples around
            int idx = rng() * activeSamples.size();
            vec3 xi = activeSamples[idx];
            bool found = false;

            // try to generate samples far away enough k times
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

            // if the area around this sample is saturated, remove this sample
            if (!found) {
                activeSamples.erase(activeSamples.begin() + idx);
#if DEBUG
                std::cout << "Removed sample. Active sample count: " << activeSamples.size() << std::endl;
#endif
            }
        }
        return gridData;
    }

    /**
     * Pairwise function to see if each sample is far away enough from all the
     * others.
     * @param allPoints
     * @return whether point set is valid poisson distribution.
     */
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

    /**
     * Create a list of points that are closer to the boundaries of the unit
     * square than the poisson disk radius.
     */
    void generateMasterTile() {
        // generate an initial poisson distribution
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

    /**
     * Generate a sample near xi using random spherical coordinates.
     * @param xi sample to generate a sample near.
     * @param r minimum distance between each pair of samples.
     * @return a new sample near xi.
     */
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

    /**
     * Check if sample's coordinates are in domain [0, 1).
     * @param sample point to test.
     * @return whether sample is in unit cube from [0, 1).
     */
    bool inGrid(vec3 sample) {
        for (int i = 0; i < 3; i++) {
            if (sample[i] < 0 || sample[i] > 1) return false;
        }
        return true;
    }

    /**
     * Check if the sample is far away enough from neighboring samples in the
     * grid.
     * @param data SamplerGrid containing all current samples.
     * @param sample candidate sample to add to the grid.
     * @param r distance to maintain between samples.
     * @return whether the candidate sample is okay to add to the grid.
     */
 	bool isFarEnough(SamplerGrid<T> data, vec3 sample, T r) {
 		vec3 gridPos = sample * inverseCellWidth;
 		vec3 gridIdx = floor(gridPos);

        // iterate over the 3x3x3 grid of neighbor cells around this sample
        int search = 1;
 		for (int z = -search; z <= search; z++) {
			for (int y = -search; y <= search; y++) {
				for (int x = -search; x <= search; x++) {
					vec3 offset(x, y, z);
					vec3 neighborIdx = offset + gridIdx;
                    int i = neighborIdx[0];
                    int j = neighborIdx[1];
                    int k = neighborIdx[2];
					if (i < 0 || i >= gridResolution ||
						j < 0 || j >= gridResolution ||
						k < 0 || k >= gridResolution)
						continue;

                    // get the list of points in the neighbor cell
                    PointList cellPts = data(i, j, k);

                    // test each point in the cell to see if its far enough
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

 	/**
 	 * Write a list of samples to a bgeo file.
 	 * @param data list of samples to write.
 	 * @param filename destination of output file.
 	 */
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

    // outer samples to use as template to generate tiles from
    std::vector<vec3> masterTile;

    // list of poisson distributions to tile
    std::vector<PointList> tileSet;

    // number of cells in each dimension of the grid
	int gridResolution;

    // 1 / cellWidth
	T inverseCellWidth;

    // width of the cells in the grid data structure
    T cellWidth;

    // dimension, poisson disk radius, and number of tries to generate successful samples
    T n, r, k;

    // number of total samples generated in this sampler
    int numSamples = 0;

    // number of poisson disk distributions to tile
    int numTiles;
};