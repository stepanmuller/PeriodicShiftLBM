PeriodicShiftLBM is a LBM solver based on the TNL library: https://tnl-project.org/

Compile using tnlcxx: https://gitlab.com/tnl-project/tnlcxx

tnlcxx --release main.cu

Features
- D3Q27 set
- periodic shift (PS) streaming
- cummulant collision taken from TNL-LBM: https://gitlab.com/tnl-project/tnl-lbm
- LES Smagorinsky model taken from: https://github.com/stloufra/LB/tree/thesis
- fully local moment-based boundary conditions
- block grid refinement
- scalar transport (in progress)

Read attached PDF for more about periodic shift, boundary conditions and grid refinement.

Performance: 2.0 GLUPS on 4070 Ti 12GB (100M cells)
