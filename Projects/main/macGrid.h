#pragma once
#include "globals.h"
#include "particle.h"
#include "gridData.h"
template <typename T>
class MACGrid {
public:
	// Constructors 
	MACGrid(int i, int j, int k) :
        particles(),
        nodeForces(i, j, k),
        nodeVels(i, j, k),
        nodeMasses(i, j, k),
        nodeMomentums(i, j, k),
        dims(3)
    {
        nodeForces.reset(vec3::Zero());
        nodeVels.reset(vec3::Zero());
        nodeMasses.reset(0);
        nodeMomentums.reset(vec3::Zero());
        dims << i, j, k;
    }

	void particleToGrid() {
        nodeVels.reset(vec3::Zero());
        nodeMasses.reset(0);
        nodeMomentums.reset(vec3::Zero());

		for (const auto &p : particles) {
            vec3 xp = p.pos;
            if (!nodeVels.inRange(xp)) { continue; }
			vec3i gridIdx = nodeVels.cellOf(xp);
            int offset = 2;

            for (int i = -offset; i <= offset; i++) {
                for (int j = -offset; j <= offset; j++) {
                    for (int k = -offset; k <= offset; k++) {
                        vec3i nIdx = gridIdx + vec3i(i, j, k);
                        if (outOfBounds(nIdx)) { continue; }
                        vec3 xi = nIdx.cast<T>() * CELL_SIZE;
                        T w = getWeight(xp, xi);
                        // transfer the masses
                        T m = w * p.mass;
                        vec3 _p = w * p.mass * (p.vel + p.B * D_INV * (xi - xp));
                        nodeMomentums(nIdx) += _p;
                        nodeMasses(nIdx) += m;
                    }
                }
            }
		}
	}


    bool outOfBounds(vec3i idx) {
        for (int m = 0; m < 3; m++) {
            if (idx[m] >= dims[m] || idx[m] < 0) {
                return true;
            }
        }
        return false;
    }
    void computeGridVelocities() {
        for (int i = 0; i < dims[0]; i++) {
            for (int j = 0; j < dims[1]; j++) {
                for (int k = 0; k < dims[2]; k++) {
                    if (nodeMasses(i, j, k) == 0) {
                        nodeVels(i, j, k) = vec3::Zero();
                    } else {
                        nodeVels(i, j, k) = nodeMomentums(i, j, k) / nodeMasses(i, j, k);
                    }
                }
            }
        }
    }

    void computeGridForces() {
        nodeForces.reset(vec3::Zero());

        for (const auto &p : particles) {
            vec3 xp = p.pos;
            if (!nodeVels.inRange(xp)) { continue; }
            vec3i gridIdx = nodeVels.cellOf(xp);
            int offset = 2;

            for (int i = -offset; i <= offset; i++) {
                for (int j = -offset; j <= offset; j++) {
                    for (int k = -offset; k <= offset; k++) {
                        vec3i nIdx = gridIdx + vec3i(i, j, k);
                        if (outOfBounds(nIdx)) continue;
                        vec3 xi = nIdx.cast<T>() * CELL_SIZE;
                        vec3 wd = getWeightGradient(xp, xi);
#ifdef STRESS
                        nodeForces(nIdx) += -V0 * NeoHookeanDeltaPsi(p.F) * p.F.transpose() * wd;
#endif
                    }
                }
            }
        }
    }


    void applyGridForces() {
        for (int i = 0; i < dims[0]; i++) {
            for (int j = 0; j < dims[1]; j++) {
                for (int k = 0; k < dims[2]; k++) {

                    nodeVels(i, j, k) += dt * nodeMasses(i, j, k) * g;
                    if (nodeMasses(i, j, k) != 0) {
                        nodeVels(i, j, k) +=  dt * nodeForces(i, j, k) / nodeMasses(i, j, k);
                    }

                }
            }
        }
    }

    void gridToParticle() {
        for (auto &p : particles) {
            vec3 xp = p.pos;
            if (!nodeVels.inRange(xp)) { continue; }
            vec3i gridIdx = nodeVels.cellOf(xp);
            int offset = 2;

            p.reset();
            mat3 newF = mat3::Zero();
            for (int i = -offset; i <= offset; i++) {
                for (int j = -offset; j <= offset; j++) {
                    for (int k = -offset; k <= offset; k++) {
                        vec3i nIdx = gridIdx + vec3i(i, j, k);
                        if (outOfBounds(nIdx)) continue;
                        vec3 xi = nIdx.cast<T>() * CELL_SIZE;
                        T w = getWeight(xp, xi);
                        vec3 v = nodeVels(nIdx);
                        newF += v * getWeightGradient(xp, xi).transpose();
                        p.vel += w * v;
#ifdef APIC
                        p.B += w * v * (xi - xp).transpose();
#else
                        p.B = mat3::Zero();
#endif
                    }
                }
            }

            p.F *= (mat3::Identity() + dt * newF);
        }
    }

    void advectParticles() {
        for (auto &p : particles) {
            p.pos += p.vel * dt;
        }
    }

    void writeFrame() {
        Partio::ParticlesDataMutable *parts = Partio::create();
        Partio::ParticleAttribute posH;
        std::string filename = "frames/frame_" + std::to_string(numFrames + 1) + ".bgeo";

        posH = parts->addAttribute("position", Partio::VECTOR, 3);

        for (unsigned int i = 0; i < particles.size(); ++i) {
            int idx = parts->addParticle();
            float *p = parts->dataWrite<float>(posH, idx);
            for (int k = 0; k < 3; ++k) {
                p[k] = particles[i].pos[k];
            }
        }

        Partio::write(filename.c_str(), *parts);
        parts->release();
        numFrames++;
    }

    void stepSimulation() {
        particleToGrid();
        computeGridVelocities();
#if DEBUG
        assert(isMomentumConserved());
#endif
        computeGridForces();
        applyGridForces();
        gridToParticle();
#if DEBUG
        assert(isMomentumConserved());
#endif
        advectParticles();
    }

    mat3 NeoHookeanDeltaPsi(mat3 F) {
        Eigen::JacobiSVD<mat3> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        mat3 U = svd.matrixU();
        mat3 V = svd.matrixV();
        if (U.determinant() < 0) { U.col(2) *= -1; }
        if (V.determinant() < 0) { V.col(2) *= -1; }

        mat3 R = U * V.transpose();

        double mu = k / (2 * (1.0 + nu));
        double lambda = (k * nu) / ((1.0 + nu) * (1.0 - (2.0 * nu)));

        mat3 JFinvT;
        JFinvT(0, 0) = F(1, 1) * F(2, 2) - F(1, 2) * F(2, 1);
        JFinvT(0, 1) = F(1, 2) * F(2, 0) - F(1, 0) * F(2, 2);
        JFinvT(0, 2) = F(1, 0) * F(2, 1) - F(1, 1) * F(2, 0);
        JFinvT(1, 0) = F(0, 2) * F(2, 1) - F(0, 1) * F(2, 2);
        JFinvT(1, 1) = F(0, 0) * F(2, 2) - F(0, 2) * F(2, 0);
        JFinvT(1, 2) = F(0, 1) * F(2, 0) - F(0, 0) * F(2, 1);
        JFinvT(2, 0) = F(0, 1) * F(1, 2) - F(0, 2) * F(1, 1);
        JFinvT(2, 1) = F(0, 2) * F(1, 0) - F(0, 0) * F(1, 2);
        JFinvT(2, 2) = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);
        double J = F.determinant();
        mat3 P = 2 * mu * (F - R) + lambda * (J - 1.f) * JFinvT;

        return P;
    }

    T getWeight(vec3 xp, vec3 xi) {
        vec3 diff = INV_CELL_SIZE * (xp - xi);
        return N(diff[0]) * N(diff[1]) * N(diff[2]);
    }

    vec3 getWeightGradient(vec3 xp, vec3 xi) {
        vec3 delta = INV_CELL_SIZE * (xp - xi);
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

    bool isMomentumConserved() {
        vec3 particleP = vec3::Zero();
        vec3 nodeP = vec3::Zero();
        for (const auto &p : particles) {
            particleP += p.mass * p.vel;
        }

        for (int i = 0; i < nodeVels.length; i++) {
            //nodeP += nodeMasses.mData[i] * nodeVels.mData[i];
            nodeP += nodeMomentums.mData[i];
        }

        bool isEqual = true;
        for (int i = 0; i < 3; i++) {
            isEqual = isEqual && std::abs(nodeP[i] - particleP[i]) < 0.00001;
            if (!isEqual) {
                printf("Particle Momentum: %f, %f, %f\n", particleP[0], particleP[1], particleP[2]);
                printf("Node Momentum: %f, %f, %f\n", nodeP[0], nodeP[1], nodeP[2]);
                return isEqual;
            }
        }
        return isEqual;
    }

	std::vector<Particle<T>> particles;
    GridData<vec3> nodeForces;
    GridData<vec3> nodeVels;
    GridData<T> nodeMasses;
    GridData<vec3> nodeMomentums;
	vec3 dims;
    int numFrames = 0;
};
