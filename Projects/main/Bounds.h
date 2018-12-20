#pragma once
#include "globals.h"

class Bounds;
struct Ray;
Bounds Union(const Bounds& b1, const Bounds& b2);
Bounds Union(const Bounds& b1, const vec3& p);
float BoundsIntersect(Bounds &b, Ray &r);

class Bounds
{
public:
    Bounds()
        : min(), max()
    {
        min << INFINITY, INFINITY, INFINITY;
        max << -INFINITY, -INFINITY, -INFINITY;
    }
    Bounds(const vec3& min, const vec3& max)
        : min(min), max(max)
    {
    }

    Bounds(const vec3& p)
        : min(p), max(p)
    {}

    Bounds(const std::vector<vec3> &pts) {
        Bounds b(pts[0]);
        for (const auto &p : pts) {
            b = Union(b, p);
        }
    }

    // Returns a vector representing the diagonal of the box
    vec3 Diagonal() const { return max - min; }

    // Returns the index of the axis with the largest length
    int MaximumExtent() const
    {
        vec3 d = Diagonal();
        if (d[0] > d[1] && d[0] > d[2])
            return 0;
        else if (d[1] > d[2])
            return 1;
        else
            return 2;
    }

    // Returns the position of a point *relative*
    // to the min and max corners of the box.
    // This ranges from (0, 0, 0) to (1, 1, 1)
    // where these two extremes represent the
    // values the result would have at the min and
    // max corners of the box
    vec3 Offset(const vec3 &p) const
    {
        vec3 o = p - min;
        if (max[0] > min[0]) o[0] /= max[0] - min[0];
        if (max[1] > min[1]) o[1] /= max[1] - min[1];
        if (max[2] > min[2]) o[2] /= max[2] - min[2];
        return o;
    }



    // Returns the surface area of this bounding box
    float SurfaceArea() const;


    //DELETEME
    inline const vec3& operator[](int i) const {
        return (i == 0) ? min : max;
    }

    //DELETEME
    inline vec3 operator[](int i) {
        return (i == 0) ? min : max;
    }


    inline bool Inside(const vec3 &p) const {
        return (p[0] >= min[0] && p[0] < max[0] &&
               p[1] >= min[1] && p[0] < max[1] &&
               p[2] >= min[2] && p[2] < max[2]);
    }

    vec3 min, max;
};
