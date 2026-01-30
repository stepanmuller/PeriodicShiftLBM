PeriodicShiftLBM is a LBM solver based on the TNL library: https://tnl-project.org/

Compile using tnlcxx: https://gitlab.com/tnl-project/tnlcxx

tnlcxx --release main.cu

Features
- periodic shift (PS) streaming
- cummulant collision taken from TNL-LBM: https://gitlab.com/tnl-project/tnl-lbm
	- just slightly modified: high order moments all get relaxed to zero so their calculation is removed
- LES taken from: https://github.com/stloufra/LB/tree/thesis
- fully local moment-based boundary conditions

Read attached PDF for more about periodic shift and BC implementation.

Performance: 2.0 GLUPS on 4070 Ti 12GB (100M cells)
