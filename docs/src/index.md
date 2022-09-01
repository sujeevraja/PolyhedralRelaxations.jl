# PolyhedralRelaxations.jl Documentation

```@meta
CurrentModule = PolyhedralRelaxations
```

## Overview

PolyhedralRelaxations.jl is a Julia/JuMP package for constructing Polyhedral Relaxations for three classes of functions (i) graphs of bounded, differentiable, univariate functions (ii) a single multilinear term where each variable involved in the multilinear term has a bounded domain and (iii) a bilinear term (special case of multilinear term). Additionally, the package also provides 4 different commonly used variable domain partitioning schemes that can be used to partition the variable domains as and when required. 

For the function class (i) it can be used to construct a sequence of both MILP and LP relaxation that converge to the graph of the univariate function and its convex hull, respectively. It takes in the function with its domain and the JuMP variables and the JuMP model and populates the model with the relaxation. It returns variable and constraint references for the formulated relaxation. 

For multilinear terms, it takes the JuMP variables and the JuMP model with the variable domains, partitions on one or more of the variable domains on each variable involved in the multilinear term and populates the model with the relaxation (if the variable domains are not partitioned, the relaxation is an LP relaxation and if not, the relaxation is an SOS2-based MILP relaxation). Similar to the univariate function relaxations, bilinear relaxation function returns the variable and constraint references for the formulated relaxation; in addition it also returns access to some additional information that were used to create the formulation. This information will in turn be useful to add strengthening constraints that link different multilinear terms that share variables between them (see later sections for description of these linking constraints). 

For bilinear term which is a special case of the bilinear term the McCormick relaxation is implemented for the LP when the variable domains are not partitioned. When exactly one of the variable domain is partitioned, an incremental formulation is implemented. This incremental formulation will work only when exactly one of the variable domains is partitioned. When both the variable domains in the bilinear term are partitioned, the multilinear relaxation can be used. 

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

## Quick Introduction to the Theory

### Geometry of relaxations for univariate functions
This quick start guide will introduce the main concepts of PolyhedralRelaxations. The package constructs MILP and LP relaxations for ``y = f(x)``, given the function and discretization points on the domain of ``f(x)``. For instance, given ``y=x^3`` and a partition of the domain ``[-1.5, 0, 2]``, it constructs an MILP relaxation of ``y = x^3`` as the disjunction of the triangles shown in the figure below:

![example](assets/example.png)

Notice from the above example that, for the above disjunction of triangles to be a valid relaxation for ``y = f(x)``, the partition points must necessarily include the inflection points of the univariate function. Also, in the above example, the domain of the function ``y = x^3`` is given by ``[-1.5, 2]``. The following table provides a list of inflection points for various well-used univariate functions. 

| Function | Inflection points |
| --- | --- |
| ``x^n`` (``n`` is odd)  | ``\{0\}`` |
| ``x^n`` (``n`` is even) | ``\phi`` |
| ``\sin x`` | ``n\pi`` ``(n \in \mathbb Z)`` | 
| ``\cos x`` | ``n\frac{\pi}{2}`` ``(n`` is odd ``)``| 
| ``x\|x\|`` | ``\{0\}``| 

If the domain of the provided function contains an inflection point, then the requirement is that the inflection points should be provided as a discretization point; failing which the MILP relaxation need not be valid for the given function. 

### Type of formulation used
The MILP formulation, obtained as the disjunction of the triangles, is formulated using the incremental method (using the so-called ``\delta`` variables) and the convex hull of the triangles (the LP) is formulated using the standard ``\lambda``-method (V-representation). For details on the formulation and theoretical properties, interested readers are referred to the following paper:

* K. Sundar, S. Sanjeevi, and H, Nagarajan (2020). Sequence of Polyhedral Relaxations for Nonlinear Univariate Functions. ([arxiv link](https://arxiv.org/abs/2005.13445))

### Geometry of relaxation for a bilinear term
Let us consider the bilinear term ``z=x_1 x_2``. Suppose ``[-1, 1]`` are the bounds on ``x_1`` and ``x_2``. Then the convex hull of ``z=x_1 x_2`` is given by the McCormick relaxation (a tetrahedron). The bilinear relaxation is then given by the four facets of the tetrahedron shown in the figure below:

![McCormick relaxation](assets/mccormick.png)

Now suppose that for the variable ``x_2``, we have partitioned the domain  has more partitions specified into ``[-1.0, -0.25, 0.25, 1.0]``, then the bilinear relaxation can be used to formulate an MILP relaxation that is given by the disjunction of three tetrahedrons shown in the figure below:

![MILP relaxation](assets/incremental-mccormick.png)

Similar to the univariate MILP relaxations, an incremental method is used to formulate this disjunction of tetrahedrons. To the best of our knowledge, this incremental formulation is not found any where in the literature and is assumed to be common knowledge. 

Currently the package does not support more than one partition on both the variables ``x_1`` and ``x_2``. Nevertheless, this support will be added in the future versions. 