[![Continuous Integration (Unit Tests)][ci-unit-img]][ci-unit-url]  [![Documentation][docs-img]][docs-url]  [![Code Coverage][codecov-img]][codecov-url]    [![Commits][commits-img]][commits-url]                                          

[docs-img]: https://github.com/sujeevraja/PolyhedralRelaxations.jl//workflows/Documentation/badge.svg "Documentation"
[docs-url]: https://kaarthiksundar.github.io/GasTranSim.jl/dev/
[ci-unit-img]: https://github.com/sujeevraja/PolyhedralRelaxations.jl/actions/workflows/ci.yml/badge.svg?branch=master "Continuous Integration (Unit Tests)"
[ci-unit-url]: https://github.com/sujeevraja/PolyhedralRelaxations.jl/actions/workflows/ci.yml
[codecov-img]: https://codecov.io/gh/sujeevraja/PolyhedralRelaxations.jl/branch/master/graph/badge.svg "Code Coverage"
[codecov-url]: https://codecov.io/gh/sujeevraja/PolyhedralRelaxations.jl/branch/master
[commits-img]: https://img.shields.io/github/commits-since/sujeevraja/PolyhedralRelaxations.jl/v0.3.5.svg "Commits since tagged version"
[commits-url]: https://github.com/sujeevraja/PolyhedralRelaxations.jl/commits/master

# PolyhedralRelaxations.jl
PolyhedralRelaxations.jl is a Julia package to construct mixed-integer linear programming and linear programming (MILP and LP) relaxations for (i) univariate, continuous, and differentiable functions whose domain is also bounded (ii) multilinear terms involving variables with bounded domain and with domain partitions on one or more variables. For bilinear terms (which is a special case of multilinear terms), we implement (using a separate API) the well known McCormick relaxation when no variable domain partition on both variables are provided and when variable domain partitions on exactly one of the variables is provided, we implement an MILP relaxation using an incremental formulation. This incremental formulation is not implemented when partitions on both variables are provided.  

## Usage

- Within the Julia REPL, run `using Pkg; Pkg.add("PolyhedralRelaxations")`
  
Check the "examples" folder for some examples how to use this package to construct relaxations for NLPs/MINLPs.

## Bug reports and support

Please report any issues via the Github **[issue tracker]**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc. 

[issue tracker]: https://github.com/sujeevraja/PolyhedralRelaxations.jl/issues

## Citation
If you find the univariate relaxations of PolyhedralRelaxations.jl useful in your work, we kindly request that you cite the following paper ([arxiv link](https://arxiv.org/abs/2005.13445)): 

```bibtex
@article{SundarSanjeeviNagarajan2021,
  title={Sequence of polyhedral relaxations for nonlinear univariate functions},
  author={Sundar, Kaarthik and Sanjeevi, Sujeevraja and Nagarajan, Harsha},
  journal={Optimization and Engineering},
  volume={23},
  number={2},
  pages={877--894},
  year={2022},
  publisher={Springer},
  doi = {10.1007/s11081-021-09609-z},
}
```

If you find the multilinear relaxations implemented in PolyhedralRelaxations.jl useful in your work, we kindly request that you cite the following two papers ([arxiv link](https://arxiv.org/abs/2001.00514) and [Optimization Online Link](https://www.optimization-online.org/DB_HTML/2022/07/8974.html)): 

```bibtex 
@article{SundarNagarajanLinderothWangBent2021,
  title={Piecewise polyhedral formulations for a multilinear term},
  author={Sundar, Kaarthik and Nagarajan, Harsha and Linderoth, Jeff and Wang, Site and Bent, Russell},
  journal={Operations Research Letters},
  volume={49},
  number={1},
  pages={144--149},
  year={2021},
  publisher={Elsevier}
}

@article{KimRichardTawarmalani2022,
    title={Piecewise Polyhedral Relaxations of Multilinear Optimization},
    author={Kim, Jongeun and Richard, Jean-Philippe P. and Tawarmalani, Mohit},
    eprinttype={Optimization Online},
    date={2022}
}
```

The MILP relaxations for the bilinear term that uses an incremental formulation is not documented any where, but is a straightforward extension of the formulation proposed for the univariate MILP relaxation in https://arxiv.org/abs/2005.13445. The MILP relaxation for nonlinear, univariate functions is a disjunction of a chain of triangles. For a bilinear term, they can be thought of as a disjunction of a chain of tetrahedrons that share edges (see the [documentation](https://sujeevraja.github.io/PolyhedralRelaxations.jl/stable/) for a visualization).





