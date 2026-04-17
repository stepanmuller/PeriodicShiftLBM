PeriodicShiftLBM is a LBM solver based on the TNL library: https://tnl-project.org/

Compile using tnlcxx: https://gitlab.com/tnl-project/tnlcxx

tnlcxx --release main.cu

Features
- D3Q27 set
- cummulant collision taken from TNL-LBM: https://gitlab.com/tnl-project/tnl-lbm
- LES Smagorinsky model taken from: https://github.com/stloufra/LB/tree/thesis
- fully local moment-based boundary conditions
- periodic shift (PS) streaming for indirectly adressed grids
- block grid refinement for indirectly adressed grids
- scalar transport (in progress, only on indirectly adressed grids with no refinement)
- esotwist streaming for directly adressed grids
- automatic multi-level grid refinement around the wall for directly adressed grids

Read attached PDF for more about periodic shift, boundary conditions and grid refinement.

Performance: 2.0 GLUPS on 4070 Ti 12GB (100M cells)
