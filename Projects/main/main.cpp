#include "globals.h"
#include "sampler.h"
#include "macGrid.h"
#include "simulation.h"

T computeSignedTetraVolume(vec3 p1, vec3 p2, vec3 p3) {
    T v321 = p3[0] * p2[1] * p1[2];
    T v231 = p2[0] * p3[1] * p1[2];
    T v312 = p3[0] * p1[1] * p2[2];
    T v132 = p1[0] * p3[1] * p2[2];
    T v213 = p2[0] * p1[1] * p3[2];
    T v123 = p1[0] * p2[1] * p3[2];
    return (1.0f / 6.0f) * (-v321 + v231 + v312 - v132 - v213 + v123);
}

T computeMeshVolume(const std::string &objPath) {
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string err;
    std::string basedir = "";

    vec3 min, max;
    min << INFINITY, INFINITY, INFINITY;
    max << -INFINITY, -INFINITY, -INFINITY;
    Bounds bounds(min, max);

    if (objPath.find_last_of("/\\") != std::string::npos) {
        basedir = objPath.substr(0, objPath.find_last_of("/\\"));
    }

    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, objPath.c_str(), basedir.c_str());
    if (!err.empty()) {
        std::cerr << err << std::endl;
    }

    if (!ret) {
        std::cerr << "Skipping " << objPath << std::endl;
        return 0.0;
    }

    T volume = 0.0;
    for (size_t s = 0; s < shapes.size(); s++) {
        // Loop over faces(polygon)
        size_t index_offset = 0;
        for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            size_t fv = 3; // all objs should have only triangles
            vec3 pts[3];
            // Loop over vertices in the face.
            for (size_t v = 0; v < fv; v++) {
                // access to vertex
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                tinyobj::real_t vx = attrib.vertices[3 * idx.vertex_index + 0];
                tinyobj::real_t vy = attrib.vertices[3 * idx.vertex_index + 1];
                tinyobj::real_t vz = attrib.vertices[3 * idx.vertex_index + 2];
                vec3 pt;
                pt << vx, vy, vz;
                pts[v] = pt;
            }


            volume += computeSignedTetraVolume(pts[0], pts[1], pts[2]);
        }
    }

    return volume;
}

int main(int argc, char **argv) {
    T DENSITY = 2.0;
    T VOLUME = 1.0;

    bool loadSampledMesh = false;
    bool loadSamples = false;
    bool loadSampledSphere = true;
	Sampler<T> *s;
	PointList samples;
    if (loadSampledMesh) {
        // load sampled mesh
        samples = loadSamplesFromFile("models/bunny.bgeo");

        // compute the bounding box
        Bounds bounds;
        for (auto &p : samples) {
            bounds = Union(bounds, p);
        }

        // move points so they are all in the positive region
        for (auto &p : samples) {
            p -= bounds.min;
        }

        VOLUME = computeMeshVolume("models/bunny.obj");
    }
	else if (loadSamples) {
        samples = loadSamplesFromFile("poisson_unit_cube_1157samples.bgeo");

	} 
	else if (loadSampledSphere) {
		//samples = loadSamplesFromFile("poisson_unit_cube_1157samples.bgeo");
        // set up sampling parameters
        T n = 3;
        T r2 = 1.0 / 3.0;
        T k = 30;

        // set up psuedo random sampling
        int seed = int(n * r2 * k);
        std::srand(seed);

        // Generating samples
        s = new Sampler<T>(n, r2 / 2.0, k);
        PointList cubeSamples = s->unitCube.getAllSamples();

        for (auto &p : cubeSamples) {
            vec3 center;
            center << 0.5, 0.5, 0.5;
            if (sphereSDF(p, center, 0.5) <= 0.0) {
                samples.push_back(p);
                vec3 higherBall = p;

                // add second sample for higher ball
                higherBall[1] += (Y_CELL_COUNT * 2 * CELL_SIZE) - 2;
                samples.push_back(higherBall);
            }
        }
        std::string size = std::to_string(samples.size());
        s->saveSamples(samples, "models/sphere_" + size + "samples.bgeo");
	}
	else {
        // set up sampling parameters
        T n = 3;
        T r2 = 1.0 / 2.0;
        T k = 30;

        // set up psuedo random sampling
        int seed = int(n * r2 * k);
        std::srand(seed);

		// Generating samples
		s = new Sampler<T>(n, r2 / 2.0, k);
		samples = s->unitCube.getAllSamples();
	}

    vec3i dim(X_CELL_COUNT, Y_CELL_COUNT, Z_CELL_COUNT);

    // create the transform for the samples
    mat4 transform = mat4::Identity();

    vec3 translate, scale;

    translate << (X_CELL_COUNT / 2) * CELL_SIZE - 0.25, 0, (Z_CELL_COUNT / 2) * CELL_SIZE - 0.25;
    scale << 0.5, 0.5, 0.5;

    for (int i = 0; i < 3; i++) {
        transform(i, 3) = translate[i];
        transform(i, i) = scale[i];
    }

    // initialize particles in the sim to have volume and mass
    std::vector<Particle<T>> particles;
    for (const auto &sample : samples) {
    	Particle<T> p(sample);
		p.vol = VOLUME / samples.size();
		p.mass = DENSITY * p.vol;
        particles.push_back(p);
    }

    // create the simulation from the parameters and particles
    Simulation sim(dim[0], dim[1], dim[2], transform, particles);

    // write the grid as a poly
    sim.grid.drawGrid();

    // run the simulation
    sim.run();
	return 0;
}
