# Yuki - 3D Elastoplastic Simulator (MPM)
## Alexander Chan, Jacob Snipes, Emily Vo

## Demos
*Insert demo videos/GIFs here*

## Implementation Details

### `globals.h`
This file is included in all other files; it stores the simulation's parameters to make for easy tuning. Further, we include custom definitions and the external Eigen and PartIO libraries.

### `particle.h`
The `Particle` class contains a position vector, velocity vector, mass quantity, volume quantity, *B* matrix as specified by the Affine Particle-In-Cell (APIC) method, elastic deformation gradient matrix (*F<sub>E</sub>*), and plastic deformation gradient matrix (*F<sub>P</sub>*).

The total deformation gradient, *F*, can then easily be computed as *F* = *F<sub>E</sub>* *F<sub>P</sub>*.

### `sampler.h`
Uses Fast Poisson Disk Sampling to generate samples within a cube. Additionally implements mesh sampling using a KD Tree (see `bounds` and `KDTree` classes) by tiling the sampled cube throughout the mesh's bounding volume and subsequently removing sampled points that lie outside the mesh.

Mesh loading was done using tinyobj.

### `gridData.h`
Template class that stores any type of data at the nodes of the marker-and-cell (MAC) grid.

### `macGrid.h`
Stores particles (of class `particles`) as well as node forces, velocities, masses, and momenta (of class `gridData`) in a 3D MAC grid and handles the main MPM algorithm:

- Transfer particle quantities to the MAC grid (using APIC)
- Compute velocities of the grid nodes
- Eliminate unnecessary grid nodes (i.e. - those with no mass) to reduce computation time
- Compute grid forces (including collision handling)
- Update grid velocities using computed forces
- Update each particle's deformation gradient
- Transfer updated grid quantities back to particles
- Advect particles appropriately
- Write data to BGEO file

This class computes the plastic deformation gradient using the first Piola-Kirchoff stress specified in Stomakhin et. al. (2013) to simulate snow-like behavior.

### `simulation.cpp`
Steps through the main algorithm and predicts time to finish writing the desired number of frames.

### `main.cpp`
Reads in OBJ files, computes total mesh volume, sets initial transformation, evenly distributes mass and volume to all particles, and runs the simulation.

## Difficulties
- Behaving like a rigid body
  - Caused by updating grid velocity with `m * g * dt` (a momentum quantity) instead of `g * dt` (a velocity quantity)
  - Prevented the body from deforming upon collision

## Future Work
- GPU parallelization using CUDA
- Scene setup and rendered simulations

## Sources
- MPM Course Notes (Jiang et. al., 2016)
  - https://www.seas.upenn.edu/~cffjiang/research/mpmcourse/mpmcourse.pdf

- "A material point method for snow simulation" (Stomakhin et. al., 2013)
  - https://www.math.ucla.edu/~jteran/papers/SSCTS13.pdf

- "Fast Poisson Disk Sampling in Arbitrary Dimensions" (Bridson, 2007)
  - https://www.cct.lsu.edu/~fharhad/ganbatte/siggraph2007/CD2/content/sketches/0250.pdf
