# PolyhedralRelaxations.jl Documentation

```@meta
CurrentModule = PolyhedralRelaxations
```

## Overview

PolyhedralRelaxations.jl is a Julia/JuMP package for constructing Polyhedral Relaxations for graphs of bounded, differentiable, univariate functions. In particular, it can be used to construct a sequence of both MILP and LP relaxation that converge to the graph of the univariate function and its convex hull, respectively. It takes in the function with its domain and returns a the constraint matrix A and a vector b that correspond to the constraint set. It also provides additional information on the type of the constraints and the variables. These can in-turn be used with a package like JuMP in an optimization model. 

## Installation Guide

To use PolyhedralRelaxations, first [download and install](https://julialang.org/downloads/) Julia or open up a remote notebook at [JuliaBox](https://www.juliabox.com/) or similar services.

This version of PolyhedralRelaxations is compatible with Julia 1.0 and later.

From Julia REPL, PolyhedralRelaxations is installed by using the built-in package manager:
```julia
import Pkg
Pkg.add("PolyhedralRelaxations")
```

## Unit Tests
To run the tests in the package, run the following command within the Julia REPL after installing the package.

```julia
import Pkg
Pkg.test("PolyhedralRelaxations")
```