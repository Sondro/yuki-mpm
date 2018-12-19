#pragma once
#include <fstream>
#include "globals.h"
#include "particle.h"
#include "gridData.h"

template <typename T>
struct SVDResult {
    mat3 U, V, Sigma;
};

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

            FOR_EACH_NEIGHBOR
                vec3i nIdx = gridIdx + vec3i(i, j, k);
                if (outOfBounds(nIdx)) { continue; }
                vec3 xi = nIdx.cast<T>() * CELL_SIZE;
                T w = getWeight(xp, xi);

                // transfer the masses
                T m = w * p.mass;
                vec3 _p = w * p.mass * (p.vel + p.B * D_INV * (xi - xp));
                nodeMomentums(nIdx) += _p;
                nodeMasses(nIdx) += m;
            END_FOR_EACH_NEIGHBOR

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

            FOR_EACH_NEIGHBOR
                vec3i nIdx = gridIdx + vec3i(i, j, k);
                if (outOfBounds(nIdx)) continue;
                vec3 xi = nIdx.cast<T>() * CELL_SIZE;
                vec3 wd = getWeightGradient(xp, xi);
                nodeForces(nIdx) += -p.vol * snowModel(p.Fe, p.Fp) * p.F.transpose() * wd;
            END_FOR_EACH_NEIGHBOR
        }
    }


    void applyGridForces() {
        for (int i = 0; i < dims[0]; i++) {
            for (int j = 0; j < dims[1]; j++) {
                for (int k = 0; k < dims[2]; k++) {

                    //nodeVels(i, j, k) += dt * nodeMasses(i, j, k) * g;
					nodeVels(i, j, k) += dt * g;
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

            p.reset();
            mat3 dws = mat3::Zero();
            FOR_EACH_NEIGHBOR
                vec3i nIdx = gridIdx + vec3i(i, j, k);
                if (outOfBounds(nIdx)) continue;
                vec3 xi = nIdx.cast<T>() * CELL_SIZE;
                T w = getWeight(xp, xi);
                vec3 v = nodeVels(nIdx);
                vec3 dw = getWeightGradient(xp, xi);
                mat3 result;
                result(0,0) = v(0) * dw(0);
                result(0,1) = v(0) * dw(1);
                result(0,2) = v(0) * dw(2);
                result(1,0) = v(1) * dw(0);
                result(1,1) = v(1) * dw(1);
                result(1,2) = v(1) * dw(2);
                result(2,0) = v(2) * dw(0);
                result(2,1) = v(2) * dw(1);
                result(2,2) = v(2) * dw(2);
                dws += result;

                p.vel += w * v;
                p.B += w * v * (xi - xp).transpose();

            END_FOR_EACH_NEIGHBOR

            // update deformation gradient
			mat3 newF = p.F + dt * dws * p.F;
			p.F = newF;

            auto Fe = p.Fe;
            auto newFe = Fe + dt * dws * Fe;

            SVDResult<T> svdRes = checkedSVD(newFe);
            mat3 Ue = svdRes.U;
            mat3 Ve = svdRes.V;
            mat3 Se = svdRes.Sigma;

            for (int i = 0; i < 3; ++i) {
                Se(i,i) = std::max(1.0 - thetaC, std::min(Se(i, i), 1.0 + thetaS));
            }

            p.Fe = Ue * Se * Ve.transpose();
            p.Fp = Ve * Se.inverse() * Ue.transpose() * newF;
        }
    }


    void advectParticles() {
        for (auto &p : particles) {
            p.pos += p.vel * dt;
        }
    }

	void setBoundaryVelocities(int thickness) {
	    for (int i = 0; i < thickness; ++i) {
	        for (int j = 0; j < dims[1]; ++j) {
	            for (int k = 0; k < dims[2]; ++k) {
	                vec3 normal;
	                normal << 1, 0, 0;
	                nodeVels(i, j, k) -= nodeVels(i, j, k).dot(normal) * normal;
	            }
	        }
	    }

	    for (int i = dims[0] - thickness; i < dims[0]; ++i) {
            for (int j = 0; j < dims[1]; ++j) {
                for (int k = 0; k < dims[2]; ++k) {
                    vec3 normal;
                    normal << -1, 0, 0;
                    nodeVels(i, j, k) -= nodeVels(i, j, k).dot(normal) * normal;
                }
            }
	    }

        for (int i = 0; i < dims[0]; ++i) {
            for (int j = 0; j < thickness; ++j) {
                for (int k = 0; k < dims[2]; ++k) {
                    vec3 normal;
                    normal << 0, 1, 0;
                    nodeVels(i, j, k) -= nodeVels(i, j, k).dot(normal) * normal;
                }
            }
        }

        for (int i = 0; i < dims[0]; ++i) {
            for (int j = dims[1] - thickness; j < dims[1]; ++j) {
                for (int k = 0; k < dims[2]; ++k) {
                    vec3 normal;
                    normal << 0, -1, 0;
                    nodeVels(i, j, k) -= nodeVels(i, j, k).dot(normal) * normal;
                }
            }
        }

        for (int i = 0; i < dims[0]; ++i) {
            for (int j = 0; j < dims[1]; ++j) {
                for (int k = 0; k < thickness; ++k) {
                    vec3 normal;
                    normal << 0, 0, 1;
                    nodeVels(i, j, k) -= nodeVels(i, j, k).dot(normal) * normal;
                }
            }
        }

        for (int i = 0; i < dims[0]; ++i) {
            for (int j = 0; j < dims[1]; ++j) {
                for (int k = dims[2] - thickness; k < dims[2]; ++k) {
                    vec3 normal;
                    normal << 0, 0, -1;
                    nodeVels(i, j, k) -= nodeVels(i, j, k).dot(normal) * normal;
                }
            }
        }

        // collision with sphere
        for (int i = 0; i < dims[0]; ++i) {
            for (int j = 0; j < dims[1]; ++j) {
                for (int k = 0; k < dims[2]; ++k) {
                    // test if cell is inside a sphere of radius 1 centered at {0, X_SIZE / 2, Z_SIZE / 2}
                    // Get world position of cell
                    vec3 pos;
                    pos << i * CELL_SIZE, j * CELL_SIZE, k * CELL_SIZE;
                    vec3 center;
                    center << X_SIZE / 2, Y_SIZE / 2, Z_SIZE / 2;
                    if ((pos - center).norm() < 0.25) {
                        // normal is pos normalize
                        vec3 normal = pos - center;
                        normal.normalize();
                        T vn = normal.dot(nodeVels(i, j, k));
                        T friction = 0.7;
                        vec3 vt = nodeVels(i, j, k) - vn * normal;
                        vec3 vt_n = vt;
                        vt_n.normalize();
                        vec3 newV = vt + friction * vn * vt_n;
                        pos.normalize();
                        nodeVels(i, j, k) =  newV;
                    }
                }
            }
        }
	}

    void writeFrame() {
        Partio::ParticlesDataMutable *parts = Partio::create();
        Partio::ParticleAttribute posH;
        Partio::ParticleAttribute velH;
        std::string filename = "frames/frame_" + std::to_string(numFrames + 1) + ".bgeo";

        posH = parts->addAttribute("position", Partio::VECTOR, 3);
        velH = parts->addAttribute("v", Partio::VECTOR, 3);

        for (unsigned int i = 0; i < particles.size(); ++i) {
            int idx = parts->addParticle();
            float *p = parts->dataWrite<float>(posH, idx);
            float *v = parts->dataWrite<float>(velH, idx);
            for (int k = 0; k < 3; ++k) {
                p[k] = particles[i].pos[k];
                v[k] = particles[i].vel[k];
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
        setBoundaryVelocities(3);
        gridToParticle();
#if DEBUG
        assert(isMomentumConserved());
#endif
        advectParticles();
    }

    mat3 NeoHookeanDeltaPsi(mat3 F) {
        mat3 U, V;
        checkedSVD(F, U, V);

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

    // Incorporates plasticity and hardening effects
    mat3 snowModel(const mat3 &Fe, const mat3 &Fp) {
        SVDResult<T> svd = checkedSVD(Fe);
        mat3 U = svd.U;
        mat3 V = svd.V;
        mat3 Re = U * V.transpose();

        // Initial values of mu and lambda; nu is Poisson's ratio
        T mu0 = k / (2 * (1.0 + nu));
        T lambda0 = (k * nu) / ((1.0 + nu) * (1.0 - (2.0 * nu)));

        // Hardening parameter; typically a value between 3 and 10
        T xi = 3.0;

        T Jp = Fp.determinant();
        T mu = mu0 * std::exp(xi * (1 - Jp));
        T lambda = lambda0 * std::exp(xi * (1 - Jp));

        T FRDet = (Fe - Re).determinant();
        T Je = Fe.determinant();
        auto piola = (2 * mu * (Fe - Re)) +
                (lambda * (Je - 1.0) *
                        Je * (Fe.transpose().inverse()));
        return piola;
    }

    SVDResult<T> checkedSVD(const mat3 &M) {
        SVDResult<T> svdRes;
        Eigen::JacobiSVD<mat3> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);
        svdRes.Sigma = mat3::Identity();
        svdRes.Sigma(0, 0) = svd.singularValues()[0];
        svdRes.Sigma(1, 1) = svd.singularValues()[1];
        svdRes.Sigma(2, 2) = svd.singularValues()[2];

        svdRes.U = svd.matrixU();
        svdRes.V = svd.matrixV();

        if (svdRes.U.determinant() < 0) {
            svdRes.U.col(2) *= -1;
            svdRes.Sigma(2, 2) *= -1;
        }

        if (svdRes.V.determinant() < 0) {
            svdRes.V.col(2) *= -1;
            svdRes.Sigma(2, 2) *= -1;
        }

        if (svdRes.Sigma(1, 1) > svdRes.Sigma(0, 0)) {
            T temp = svdRes.Sigma(1, 1);
            svdRes.Sigma(1, 1) = svdRes.Sigma(0, 0);
            svdRes.Sigma(0, 0) = temp;
        }

        return svdRes;
    }

    void checkedSVD(const mat3 &F, mat3& U, mat3& V) {
        Eigen::JacobiSVD<mat3> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        U = svd.matrixU();
        V = svd.matrixV();
        if (U.determinant() < 0) { U.col(2) *= -1; }
        if (V.determinant() < 0) { V.col(2) *= -1; }
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
            return 0.5 * x * (3 * absX - 4.0);
        } else if (absX >= 1.0 && absX < 2.0) {
            T k = 2.0 - absX;
            return (-x * k * k) / (2 * absX);
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
            nodeP += nodeMasses.mData[i] * nodeVels.mData[i];
            //nodeP += nodeMomentums.mData[i];
        }

        bool isConserved = cmpVec3(nodeP, "Node Momentums", particleP, "Particle Momentums");
        return isConserved;
    }

    bool cmpVec3(vec3 &v1, const char *v1Name, vec3 &v2, const char *v2Name) {
        bool isEqual = true;
        for (int i = 0; i < 3; i++) {
            isEqual = isEqual && std::abs(v1[i] - v2[i]) < EPSILON;
            if (!isEqual) {
                printf("%s: %f, %f, %f\n", v1Name, v1[0], v1[1], v1[2]);
                printf("%s: %f, %f, %f\n", v2Name, v2[0], v2[1], v2[2]);
                return isEqual;
            }
        }
        return isEqual;
    }

    bool cmpT(T v1, const char *v1Name, T v2, const char *v2Name) {
        if (std::abs(v1 - v2) > EPSILON) {
            printf("%s: %f", v1Name, v1);
            printf("%s: %f", v2Name, v2);
            return false;
        } else {
            return true;
        }
    }

    void drawGrid() {
        // dims is number of cells along an axis
        std::ofstream polyFile("grid.poly");

        polyFile << "POINTS\n";
        for (unsigned int k = 0; k < dims[2]; ++k) {
            for (unsigned int j = 0; j < dims[1]; ++j) {
                for (unsigned int i = 0; i < dims[0]; ++i) {
                    unsigned int idx = i + j * dims[0] + k * dims[0] * dims[1];
                    polyFile << idx + 1 << ": ";
                    T x = CELL_SIZE * i;
                    T y = CELL_SIZE * j;
                    T z = CELL_SIZE * k;
                    polyFile << x << " " << y << " " << z << "\n";
                }
            }
        }

        polyFile << "POLYS\n";
        int segCount = 1;
        for (int i = 0; i < dims[0]; ++i) {
            for (int j = 0; j < dims[1]; ++j) {
                for (int k = 0; k < dims[2]; ++k) {
                    int idx = i + j * dims[0] + k * dims[0] * dims[1];
                    int idx_0 = (i - 1 == -1 ? 0 : i - 1) +
                                 j * dims[0] +
                                 k * dims[0] * dims[1];
                    int idx_1 = (i + 1 == dims[0] ? dims[0] - 1 : i + 1) +
                                 j * dims[0] +
                                 k * dims[0] * dims[1];
                    int idx_2 = i +
                                 (j - 1 == -1 ? 0 : j - 1) * dims[0] +
                                 k * dims[0] * dims[1];
                    int idx_3 = i +
                                 (j + 1 == dims[1] ? dims[1] : j + 1) * dims[0] +
                                 k * dims[0] * dims[1];
                    int idx_4 = i +
                                 j * dims[0] +
                                 (k - 1 == -1 ? 0 : k - 1) * dims[0] * dims[1];
                    int idx_5 = i +
                                 j * dims[0] +
                                 (k + 1 == dims[2] ? dims[2] : k + 1) * dims[0] * dims[1];
                    int indices[6] = { idx_0, idx_1, idx_2, idx_3, idx_4, idx_5 };
                    for (int l = 0; l < 6; ++l) {
                        polyFile << segCount << ":";
                        polyFile << " " << idx << " " << indices[l] << "\n";
                        segCount++;
                    }
                }
            }
        }

        polyFile << "END";
        polyFile.close();
    }

	std::vector<Particle<T>> particles;
    GridData<vec3> nodeForces;
    GridData<vec3> nodeVels;
    GridData<T> nodeMasses;
    GridData<vec3> nodeMomentums;
	vec3 dims;
    int numFrames = 0;
};
