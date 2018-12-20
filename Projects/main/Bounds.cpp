#include "Bounds.h"
#include "KDTree.h"

float Bounds::SurfaceArea() const
{
    vec3 d = Diagonal();
    return 2 * (d[0] * d[1] + d[0] * d[2] + d[1] * d[2]);
}

Bounds Union(const Bounds& b1, const Bounds& b2)
{
    return Bounds(vec3(std::min(b1.min[0], b2.min[0]),
                            std::min(b1.min[1], b2.min[1]),
                            std::min(b1.min[2], b2.min[2])),
                    vec3(std::max(b1.max[0], b2.max[0]),
                            std::max(b1.max[1], b2.max[1]),
                            std::max(b1.max[2], b2.max[2])));
}

Bounds Union(const Bounds& b1, const vec3& p)
{
    return Bounds(vec3(std::min(b1.min[0], p[0]),
                            std::min(b1.min[1], p[1]),
                            std::min(b1.min[2], p[2])),
                    vec3(std::max(b1.max[0], p[0]),
                            std::max(b1.max[1], p[1]),
                            std::max(b1.max[2], p[2])));
}


float BoundsIntersect(Bounds &b, Ray &r) {
    float tmin = -1e38f;
    float tmax = 1e38f;
    for (int xyz = 0; xyz < 3; ++xyz) {
        float qdxyz = r.direction[xyz];
        /*if (glm::abs(qdxyz) > 0.00001f)*/
        {
            float t1 = (b.min[xyz] - r.origin[xyz]) / qdxyz;
            float t2 = (b.max[xyz] - r.origin[xyz]) / qdxyz;
            float ta = std::min(t1, t2);
            float tb = std::max(t1, t2);
            if (ta > 0 && ta > tmin)
                tmin = ta;
            if (tb < tmax)
                tmax = tb;
        }
    }
    if (tmax >= tmin && tmax > 0) {
        return tmin;
    }
    return -1;
}
