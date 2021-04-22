[![Build Status](https://github.com//sujeevraja/PolyhedralRelaxations.jl/workflows/CI/badge.svg?branch=master)](https://github.com/sujeevraja/PolyhedralRelaxations.jl/actions?query=workflow%3ACI) 
[![codecov](https://codecov.io/gh/sujeevraja/PolyhedralRelaxations.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/sujeevraja/PolyhedralRelaxations.jl)

Latest: [![](https://github.com//sujeevraja/PolyhedralRelaxations.jl/workflows/Documentation/badge.svg)](https://sujeevraja.github.io/PolyhedralRelaxations.jl/stable/)
Dev: [![](https://github.com//sujeevraja/PolyhedralRelaxations.jl/workflows/Documentation/badge.svg)](https://sujeevraja.github.io/PolyhedralRelaxations.jl/dev/)

# PolyhedralRelaxations.jl
PolyhedralRelaxations.jl is a Julia package to construct mixed-integer linear programming and linear programming (MILP and LP) relaxations for (i) univariate, continuous, and differentiable functions whose domain is also bounded (ii) bilinear terms with partitions on one of the variables.

## Usage

- Clone the repository.
- Open a terminal in the repo folder and run `julia --project=.`.
- Hit `]` to open the project environment and run `test` to run unit tests. If
  you see an error because of missing packages, run `resolve`.

Check the "examples" folder on how to use this package.

## Bug reports and support

Please report any issues via the Github **[issue tracker]**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc. 

[issue tracker]: https://github.com/sujeevraja/PolyhedralRelaxations.jl/issues

## Citation
If you find PolyhedralRelaxations.jl useful in your work, we kindly request that you cite the following paper ([arxiv link](https://arxiv.org/abs/2005.13445)): 

If you find JuMP useful in your work, we kindly request that you cite the
following paper ([pdf](https://mlubin.github.io/pdf/jump-sirev.pdf)):

```bibtex
@article{SundarSanjeeviNagarajan2021,
  title={Sequence of polyhedral relaxations for nonlinear univariate functions},
  author={Sundar, Kaarthik and Sanjeevi, Sujeevraja and Nagarajan, Harsha},
  journal={Optimization and Engineering},
  pages={1--18},
  year={2021},
  publisher={Springer},
  doi = {10.1007/s11081-021-09609-z},
}
```

The MILP relaxations for the bilinear term is not documented any where, but is an extension of the formulation proposed for the univariate MILP relaxation in the above paper. The MILP relaxation for nonlinear, univariate functions is a disjunction of a chain of triangles. For a bilinear term, they can be thought of as a disjunction of a chain of tetrahedrons that share edges (see the [documentation](https://sujeevraja.github.io/PolyhedralRelaxations.jl/stable/) for a visualization).





