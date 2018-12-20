#pragma once

#include "globals.h"
#include "Bounds.h"
#include <unordered_set>

struct Triangle {
	Triangle() {}
	Triangle(vec3 p1, vec3 p2, vec3 p3) {
		pts[0] = p1;
		pts[1] = p2;
		pts[2] = p3;
		bounds = Bounds({p1, p2, p3});
	}
	vec3 pts[3];
	Bounds bounds;
};

struct Ray {
	vec3 origin;
	vec3 direction;
};

inline float getArea(vec3 p1, vec3 p2, vec3 p3) {
	vec3 AB = p2 - p1;
	vec3 AC = p3 - p1;
	return 0.5 * AB.cross(AC).norm();
}

inline bool triangleIntersectionTest(Triangle &tri, Ray &r) {
	vec3 e1 = tri.pts[0] - tri.pts[1];
	vec3 e2 = tri.pts[2] - tri.pts[1];
    vec3 normal = -e1.cross(e2);
    normal.normalize();
    float D = normal.dot(tri.pts[0]);
    if (std::abs(D) < 0.001) { return false; }
    float t = normal.dot(tri.pts[0] - r.origin) / normal.dot(r.direction);
    if (t < 0) { return false; }
    vec3 P = r.origin + t * r.direction;
	float S = getArea(tri.pts[0], tri.pts[1], tri.pts[2]);
	float S1 = getArea(P, tri.pts[1], tri.pts[2]) / S;
	float S2 = getArea(P, tri.pts[2], tri.pts[0]) / S;
	float S3 = getArea(P, tri.pts[0], tri.pts[1]) / S;
	return (0 <= S1 && S1 <= 1 &&
			0 <= S2 && S2 <= 1 &&
			0 <= S3 && S3 <= 1 &&
			(S1 + S2 + S3 - 1 < 0.001));
}


class KDNode {
public:
	KDNode(std::vector<Triangle *> &pts, int depth, int maxDepth = 20) {
		left = nullptr;
		right = nullptr;
		bounds = pts[0]->bounds;
		for (const auto &t : pts) {
			bounds = Union(bounds, t->bounds);
		}
		if (pts.size() <= 4 || depth > maxDepth) {
			tris = pts;
			return;
		}

		int axis = bounds.MaximumExtent();
		if (std::abs(bounds.max[axis] - bounds.min[axis]) < 0.001) {
			tris = pts;
			return;
		}
		std::vector<Triangle *> rightGeos;
		std::vector<Triangle *> leftGeos;

		std::sort(pts.begin(), pts.end(), [&](const Triangle *t1, const Triangle *t2) {
			return (t1->bounds.Diagonal() / 2.0)[axis] < (t2->bounds.Diagonal() / 2.0)[axis];
		});
		leftGeos = std::vector<Triangle *>(pts.begin(), pts.begin() + (pts.size() / 2));
		rightGeos = std::vector<Triangle *>(pts.begin() + (pts.size() / 2), pts.end());
		left = new KDNode(leftGeos, depth + 1, maxDepth);
		right = new KDNode(rightGeos, depth + 1, maxDepth);
		return;
	}

	int countIntersections(Ray &r, std::unordered_set<Triangle *> &checked) {
		int totalIntersections = 0;
		if (BoundsIntersect(bounds, r) == -1) {
			return 0;
		}
		if (left == nullptr || right == nullptr) {
			// Leaf, check and count intersection with each triangle
			for (auto &t : tris) {
                if (checked.find(t) != checked.end()) { continue; }
				if (triangleIntersectionTest(*t, r)) {
					totalIntersections++;
				}
                checked.emplace(t);
			}
			return totalIntersections;
		}
		// Check bounds intersect with this node
		// If true, then count intersections in left and right (if not nullptr)
		totalIntersections += left->countIntersections(r, checked);
		totalIntersections += right->countIntersections(r, checked);
		return totalIntersections;
	}

	Bounds bounds;
	KDNode *left;
	KDNode *right;
	std::vector<Triangle *> tris;
};