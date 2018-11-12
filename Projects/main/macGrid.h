#pragma once
#include "globals.h"
#include "particle.h"
#include "gridData.h"
template <typename T>
class MACGrid {
public:
	// Constructors 
	MACGrid(vec3i dims) :
        particles(),
        nodeForces(dims),
        nodeVels(dims),
        nodeMasses(dims),
        nodeMomentums(dims),
        dims(dims)
    {}

	void particleToGrid() {
        nodeVels.reset();
		for (const auto &p : particles) {
            vec3 xp = p.pos;
			vec3i gridIdx = nodeVels.cellOf(xp);
            int offset = 2;

            for (int i = -offset; i <= offset; i++) {
                for (int j = -offset; j <= offset; j++) {
                    for (int k = -offset; k <= offset; k++) {
                        vec3 nIdx = gridIdx + vec3(i, j, k);
                        vec3 xi = nIdx * CELL_SIZE;
                        T w = getWeight(xp, xi);
                        // transfer the masses
                        T m = w * p.mass;
                        T p = w * p.mass * (p.vel + p.B * D_INV * (xi - xp));
                        nodeMomentums(nIdx) += p;
                        nodeMasses(nIdx) += m;
                    }
                }
            }
		}
	}

    void computeGridVelocities() {
        for (int i = 0; i < dims[0]; i++) {
            for (int j = 0; j < dims[1]; j++) {
                for (int k = 0; k < dims[2]; k++) {
                    nodeVels(i, j, k) = nodeMomentums(i, j, k) / nodeMasses(i, j, k);
                }
            }
        }
    }

    void computeGridForces() {
        nodeForces.reset();
        for (const auto &p : particles) {
            vec3 xp = p.pos;
            vec3i gridIdx = nodeVels.cellOf(xp);
            int offset = 2;

            for (int i = -offset; i <= offset; i++) {
                for (int j = -offset; j <= offset; j++) {
                    for (int k = -offset; k <= offset; k++) {
                        vec3 nIdx = gridIdx + vec3(i, j, k);
                        vec3 xi = nIdx * CELL_SIZE;
                        vec3 wd = getWeightGradient(xp, xi);
                        nodeForces(nIdx) += -V0 * NeoHookeanDeltaPsi(p.F) * p.F.transpose() * wd;
                    }
                }
            }
        }
    }


    void applyGridForces() {
        for (int i = 0; i < dims[0]; i++) {
            for (int j = 0; j < dims[1]; j++) {
                for (int k = 0; k < dims[2]; k++) {
                    nodeVels(i, j, k) += dt * nodeMasses(i, j, k) * 9.8;
                    nodeVels(i, j, k) += dt * nodeForces(i, j, k) / nodeMasses(i, j, k);
                }
            }
        }
    }

    void gridToParticle() {
        for (const auto &p : particles) {
            vec3 xp = p.pos;
            vec3i gridIdx = nodeVels.cellOf(xp);
            int offset = 2;

            p.reset();
            mat3 newF(0);
            for (int i = -offset; i <= offset; i++) {
                for (int j = -offset; j <= offset; j++) {
                    for (int k = -offset; k <= offset; k++) {
                        vec3 nIdx = gridIdx + vec3(i, j, k);
                        vec3 xi = nIdx * CELL_SIZE;
                        T w = getWeight(xp, xi);
                        T v = nodeVels(nIdx);
                        newF += v * getWeightGradient(xp, xi).transpose();
                        p.vel += w * v;
                        p.B += w * v * (xi - xp).transpose();
                    }
                }
            }

            p.F *= (mat3::Identity() + dt * newF);
        }
    }

    void advectParticles() {
        for (const auto &p : particles) {
            p.pos += dt * p.vel;
        }
    }

    void writeFrame() {
        Partio::ParticlesDataMutable *parts = Partio::create();
        Partio::ParticleAttribute posH;
        std::string filename = "frames/frame_" + std::to_string(numFrames) + ".bgeo";

        posH = parts->addAttribute("position", Partio::VECTOR, 3);

        for (int i = 0; i < particles.size(); ++i) {
            int idx = parts->addParticle();
            float *p = parts->dataWrite<float>(posH, idx);
            for (int k = 0; k < 3; ++k) {
                p[k] = particles[i][k];
            }
        }

        Partio::write(filename.c_str(), *parts);
        parts->release();
        numFrames++;
    }

    void stepSimulation() {
        particleToGrid();
        computeGridVelocities();
        computeGridForces();
        applyGridForces();
        gridToParticle();
        advectParticles();
        writeFrame();
    }

    mat3 NeoHookeanDeltaPsi(mat3 F) {
        double mu = k / (2 * (1.0 + nu));
        double lambda = (k * nu) / ((1.0 + nu) * (1.0 - (2.0 * nu)));
        mat3 E = (1.0 / 2.0) * (F.transpose() * F - mat3::Identity());
        mat3 F_invTrans = F.inverse().transpose();
        float J = F.determinant();
        float logJ = std::log(J);
        mat3 P = mu * (F - mu * F_invTrans) + lambda * logJ * F_invTrans;
        return P;
    }

    T getWeight(vec3 xp, vec3 xi) {
        vec3 diff = INV_CELL_SIZE * (xp - xi);
        return N(diff[0]) * N (diff[1]) * N(diff[2]);
    }

    vec3 getWeightGradient(vec3 xp, vec3 xi) {
        T delta = INV_CELL_SIZE * (xp - xi);
        vec3 n(N(delta[0]), N(delta[1]), N(delta[2]));
        vec3 dn(NPrime(delta[0]), NPrime(delta[1]), NPrime(delta[2]));
        vec3 wd;
        wd[0] = INV_CELL_SIZE * dn[0] * n[1] * n[2];
        wd[1] = n[0] * INV_CELL_SIZE * dn[1] * n[2];
        wd[2] = n[0] * n[1] * INV_CELL_SIZE * dn[2];
        return wd;
    }

    // Cubic interpolation kernel
    T N(T x) {
        T absX = std::abs(x);
        if (absX >= 0.0 && absX < 1.0) {
            T x2 = absX * absX;
            T x3 = x2 * absX;
            return 0.5 * x3 - x2 + (2.0 / 3.0);
        } else if (absX >= 1.0 && absX < 2.0) {
            T k = 2.0 - absX;
            return (1.0 / 6.0) * k * k * k;
        } else {
            return 0.0;
        }
    }

    // Cubic interpolation kernel derivative
    T NPrime(T x) {
        T absX = std::abs(x);
        if (absX >= 0.0 && absX < 1.0) {
            T x2 = absX * absX;
            return 1.5 * x2 - 2 * absX;
        } else if (absX >= 1.0 && absX < 2.0) {
            T k = 2.0 - absX;
            return -0.5 * k * k;
        } else {
            return 0.0;
        }
    }

	std::vector<Particle<T>> particles;
    GridData<vec3> nodeForces;
    GridData<vec3> nodeVels;
    GridData<T> nodeMasses;
    GridData<vec3> nodeMomentums;
	vec3i dims;
    int numFrames = 0;
};
