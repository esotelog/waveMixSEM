# Introduction
This code uses mix CG-DG (SEM) formulations to embed  a large fracture (between DG cells) in a homogeneous background (CG cells). The effect of the fracture is included only via a special flux but not physically on the mesh. The method solves for the elastic wave equation.

# To use this code
./waveMixSEM_exe param_waveMixSEM.prm

# example folder
## Run the code
./waveMixSEM_exe example.prm

output
- example.dat
- vtk files 
- waveMixSEM_log.txt

## Model: Shear plane wave
![Simulation model](./example/figures/model_example.png)

- Transient shear traction with ricker wavelet.
- Right and left periodic boundary conditions.
- Top and bottom first order absorbing boundary conditions.

## Seismograms
1. R1: top receiver
![Top  receiver](./example/figures/example_topreceiver.png)

Top receiver showing direct and reflected waves.

2. R2: bottom receiver
![Top  receiver](./example/figures/example_bottomreceiver.png)

Bottom receiver showing trasmitted wave

## Animation

![Wave propagation simulation](./example/figures/example.gif)

[▶ Download full video (MP4)](./example/figures/example.mp4)
