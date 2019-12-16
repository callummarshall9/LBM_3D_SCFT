# LBM with SCFT
Implementation of Lattice Boltzmann in 3D using D3Q15 that also solves
the field equations for a polymer chain.

# Lattice Boltzmann Implementation

1. Implementation, initialises using the equilibrium distribution.

2. Time Step Algorithm
- Macroscopic moments are computed
- Equilibrium distribution is calculated.
- The macroscopic fields density and fluid speed is written to hard disk for visualisation or post processing.
- Colllision relaxation is performed.
- Streaming is done.
- The time step is incremented and repeat the time step algorithm.

Implemeted using The Lattice Boltzmann Method Principles and Practice. ISBN: 978-3-319-44647-9
SCFT solver implemented using "Lattice Boltzmann method for multiscale self-consistent field theory
simulations of block copolymers" by: Hsieh Chen, YongJoo Kim, and Alfredo Alexander-Katz
and "Lattice Boltzmann study of pattern formation in reaction-diffusion" by: S.G. Ayodele1, F. Varnik
and D. Raabe

Lees-Edwards boundary conditions implemented using: "Lees-Edwards boundary conditions for lattice Boltzmann suspension simulations" by: Eric Lorenz and Alfons G. Hoekstra. As well as: "Lees-Edwards boundary conditions for Lattice
Boltzmann" by: Alexander J. Wagner and Ignacio Pagonabarraga.

# To build

To build. Clone repo. Then
mkdir build \
cd build \
cmake .. \
make

# Projects used in this
- RapidJSON
- fast-cpp-csv-parser (https://github.com/ben-strasser/fast-cpp-csv-parser)

# To visualise
Run the simulation, then once it has ran. Run the visualise MATLAB script and it will then generate figure images and interactive figures in MATLAB.

# Notes

options.json must be copied to the build folder currently. \
The options.json has save_every option which generates a CSV file of the density and velocity fields for the current timestep % save_every. So a save_every of 550 will generate 0.csv,50.csv,100.csv and so forth. Grid size is the grid of the field to visualise size. Velocity set can either be D3Q15 or D3Q27. c_s can be adjusted as different velocity sets may have alternate speed of sounds.

# Cloning instructions

git clone https://github.com/callummarshall9/LBM_3D_SCFT.git \
git submodule update --init --recursive


