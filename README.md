Build status: [![Build Status](https://travis-ci.org/sujeevraja/PolyhedralRelaxations.jl.svg?branch=master)](https://travis-ci.org/sujeevraja/PolyhedralRelaxations.jl) 

Code coverage: [![codecov](https://codecov.io/gh/sujeevraja/PolyhedralRelaxations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/sujeevraja/PolyhedralRelaxations.jl)
[![Coverage Status](https://coveralls.io/repos/github/sujeevraja/PolyhedralRelaxations.jl/badge.svg?branch=master)](https://coveralls.io/github/sujeevraja/PolyhedralRelaxations.jl?branch=master)

Documentation: [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://sujeevraja.github.io/PolyhedralRelaxations.jl/stable/)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://sujeevraja.github.io/PolyhedralRelaxations.jl/latest/)

# PolyhedralRelaxations.jl
PolyhedralRelaxations.jl is a Julia package to construct mixed-integer linear programming and linear programming (MILP and LP) relaxations for univariate, continuous, and differentiable functions whose domain is also bounded. The package constructs the relaxations and returns a sparse matrix A and column vector b with additional details on what each variable represents and whether the variable is continuous or binary. 

## Bug reports and support

Please report any issues via the Github **[issue tracker]**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc. 

[issue tracker]: https://github.com/sujeevraja/PolyhedralRelaxations.jl/issues

## Citation
If you find PolyhedralRelaxations.jl useful in your work, we kindly request that you cite the following paper: 

K. Sundar, S. Sanjeevi, and H, Nagarajan (2020). Sequence of Polyhedral Relaxations for Nonlinear Univariate Functions. ([arxiv link](https://arxiv.org/abs/2005.13445))





