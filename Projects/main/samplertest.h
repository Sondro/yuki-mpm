#pragma once

#include "globals.h"
#include "sampler.h"

template <typename T>
class SamplerTest {
public:
    SamplerTest(Sampler<T> *s) : sampler(s) {}

    bool testUnitPoissonCube() {
        SamplerGrid<T> result = sampler->generatePoissonDistr();
        PointList samples = result.getAllSamples();

        sampler->saveSamples(samples, "fixedPoisson.bgeo");
        return sampler->validPointSet(samples);
    }

    bool testTileCubeDistribution() {
        std::cout << "Test tiling in x and y directions..." << std::endl;
        Partio::ParticlesDataMutable *parts = Partio::create();
        Partio::ParticleAttribute posH;

        posH = parts->addAttribute("position", Partio::VECTOR, 3);
        std::string filename = "multipleTiles.bgeo";
        std::vector<vec3> allPoints;
        for (int n = 0; n < sampler->numTiles; n++) {
            std::vector<vec3> data = sampler->tileSet[n];
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

        return sampler->validPointSet(allPoints);
    }
    bool testTileRowDistribution() {
        std::cout << "Test tile row distribution..." << std::endl;
        Partio::ParticlesDataMutable *parts = Partio::create();
        Partio::ParticleAttribute posH;

        posH = parts->addAttribute("position", Partio::VECTOR, 3);
        std::string filename = "multipleTiles.bgeo";
        std::vector<vec3> allPoints;
        for (int n = 0; n < sampler->numTiles; n++) {
            std::vector<vec3> data = sampler->tileSet[n];
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

        return sampler->validPointSet(allPoints);
    }

    Sampler<T> *sampler;
};