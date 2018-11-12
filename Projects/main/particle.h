#pragma once
#include "globals.h"

template <typename T>
class Particle {
public:
    Particle(vec3 pos = vec3(0),
             vec3 vel = vec3(0),
             T mass = 1.0,
             mat3 B = mat3::Zero(),
             mat3 F = mat3::Identity()) :
             pos(pos), vel(vel), mass(mass), B(B), F(F) {}
	vec3 pos;
	vec3 vel;
	T mass;
	mat3 B;
    mat3 F; // deformation gradient
    void reset() {
        this->vel = vec3(0);
        this->B = mat3(0);
    }
};
