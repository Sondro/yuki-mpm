#pragma once

#include "globals.h"

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
        if (idx < 0 || static_cast<unsigned int>(idx) >= mData.size()) return false;
        (*this)[idx].push_back(pt);
        return true;
    }

    /**
     * Return a list of samples contained in this grid data structure.
     * @return lists of samples.
     */
    PointList getAllSamples() const {
    	if (!isAllSamplesGen) {
    		allSamples.clear();
			for (int i = 0; i < gridResolution; i++) {
				for (int j = 0; j < gridResolution; j++) {
					for (int k = 0; k < gridResolution; k++) {
						PointList pointList = (*this)(i, j, k);
						for (vec3 p : pointList) {
							allSamples.push_back(p);
						}
					}
				}
			}
    	}
    	return allSamples;
    }

    int gridResolution; // number of cells for each side.
    T cellWidth; // width of the cells in the grid.
    T inverseCellWidth; // 1 / cellWidth.

private:
    FlatGrid mData; // array of point lists. there is a list of points per cell.
    mutable bool isAllSamplesGen;
    mutable PointList allSamples;
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
    Sampler(T n, T r, T k) :
        n(n), r(r), k(k),
        gridResolution(1.0 / (2.0 * r)),
        inverseCellWidth(1.0 / (2.0 * r)),
        cellWidth(2.0 * r),
        numSamples(0),
        unitCube(generatePoissonDistr()) {}

	Sampler(T n, T r, T k, const std::string &filename) :
		n(n), r(r), k(k),
		gridResolution(1.0 / (2.0 * r)),
		inverseCellWidth(1.0 / (2.0 * r)),
		cellWidth(2.0 * r),
		numSamples(0),
		unitCube(loadSamples(filename)) {}

    /**
     * Generate poisson samples using Fast Poisson Disk Sampling. Grid is
     * initialized with a single sample.
     * @return a sampler grid data structure filled with poisson disk samples.
     */
 	SamplerGrid<T> generatePoissonDistr() {
        // initialize 3D grid
        SamplerGrid<T> gridData(gridResolution, cellWidth);

        // add first sample to the grid
        vec3 x0(rng(), rng(), rng());
        gridData.addSample(x0);

        // make first sample active
        std::vector<vec3> activeSamples;
        activeSamples.push_back(x0);

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
                        printf("Added sample: %f, %f, %f\n", sample[0], sample[1], sample[2]);
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
    bool validPointSet(PointList &allPoints) {
        bool testPassed = true;
        int numClosePairs = 0;
        for (int i = 0; i < allPoints.size(); ++i) {
            for (int j = i + 1; j < allPoints.size(); ++j) {
                vec3 diff = allPoints[j] - allPoints[i];
                if (diff.norm() - r < -0.001) {
                    std::cout << "Points too close!" << std::endl;
                    std::cout << "Distance: " << diff.norm() << std::endl;
                    std::cout << "Poisson Radius: " << r << std::endl;
                    printf("Pos 1: %f, %f, %f\n", allPoints[j][0], allPoints[j][1], allPoints[j][2]);
                    printf("Pos 2: %f, %f, %f\n", allPoints[i][0], allPoints[i][1], allPoints[i][2]);
                    printf("i: %d, j: %d\n", i, j);
                    testPassed = false;
                    numClosePairs += 1;
                }
            }
        }

        if (testPassed) std::cout << "Points meet poisson disk radius constraint." << std::endl;
        else {
            std::cout << "Number of close pairs: " << numClosePairs << std::endl;
        }
        return testPassed;
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
                    vec3 wrappedSample = sample;
					vec3 offset(x, y, z);
					vec3 neighborIdx = offset + gridIdx;

                    // if the neighbor exceeds the bounds of the unit cube
                    // wrap the neighbor index and the point we're testing
                    for (int n = 0; n < 3; n++) {
                        if (int(neighborIdx[n]) >= gridResolution) {
                            neighborIdx[n] = 0;
                            wrappedSample[n] -= 1;
                        } else if (int(neighborIdx[n]) < 0) {
                            neighborIdx[n] = gridResolution - 1;
                            wrappedSample[n] += 1;
                        }
                    }
                    int i = neighborIdx[0];
                    int j = neighborIdx[1];
                    int k = neighborIdx[2];

                    // get the list of points in the neighbor cell
                    PointList cellPts = data(i, j, k);

                    // test each point in the cell to see if its far enough
                    if (!cellPts.empty()) {
                        for (vec3 &neighborSample : cellPts) {
                            vec3 v = neighborSample - wrappedSample;
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

 	SamplerGrid<T> loadSamples(const std::string &filename) {
 		std::cout << "Loading samples from " << filename << "\n";
		SamplerGrid<T> samples(gridResolution, cellWidth);
 		Partio::ParticlesDataMutable *parts = Partio::read(filename.c_str());
 		Partio::ParticleAttribute posAttr;
 		if (!parts->attributeInfo("position", posAttr)) {
 			std::cerr << "Invalid sample set\n";
 		}

 		for (int i = 0; i < parts->numParticles(); ++i) {
 			const float *pos = parts->data<float>(posAttr, i);
 			vec3 newPt;
 			newPt << pos[0], pos[1], pos[2];
 			samples.addSample(newPt);
 		}
 		std::cout << "Loaded " << samples.getAllSamples().size() << " samples\n";
 		return samples;
 	}

    // dimension, poisson disk radius, and number of tries to generate successful samples
    T n, r, k;

    // number of cells in each dimension of the grid
	int gridResolution;

    // 1 / cellWidth
	T inverseCellWidth;

    // width of the cells in the grid data structure
    T cellWidth;

    // number of total samples generated in this sampler
    int numSamples = 0;

    // grid full of points in domain [0, 1)
    SamplerGrid<T> unitCube;
};