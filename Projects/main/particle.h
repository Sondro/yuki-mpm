#pragma once
#include "globals.h"

template <typename T>
class Particle {
public:
    Particle(vec3 pos = vec3::Zero(),
             vec3 vel = vec3::Zero(),
             T mass = 1.0,
             mat3 B = mat3::Zero(),
             mat3 F = mat3::Identity(),
             mat3 Fe = mat3::Identity(),
             mat3 Fp = mat3::Identity()) :
             pos(pos), vel(vel), mass(mass), B(B), F(F), Fe(Fe), Fp(Fp) {}
	vec3 pos;
	vec3 vel;
	T mass;
	mat3 B;
    mat3 F; // deformation gradient; F = Fe * Fp
    mat3 Fe; // Elastic deformation
    mat3 Fp; // Plastic deformation
    T vol;
    void reset() {
        this->vel = vec3::Zero();
        this->B = mat3::Zero();
    }
};
