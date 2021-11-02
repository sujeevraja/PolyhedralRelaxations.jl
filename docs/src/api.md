API
=================


## Basic Usage
Here is an example illustrating the basic API for univariate functions.

```julia 
using PolyhderalRelaxations, JuMP
m = Model() 
@variable(m, -1.0 <= x <= 1.0)
@variable(m, y)
f = x -> x^3
construct_univariate_relaxation!(m, f, x, y, [-1.0, 0.0, 1.0], true)
```

For bilinear functions, the basic API is as follows:

```julia 
using PolyhderalRelaxations, JuMP
m = Model() 
@variable(m, -1.0 <= x <= 1.0)
@variable(m, -10.0 <= y <= 12.0)
@variable(m, z)
construct_bilinear_relaxation!(m, x, y, z, [-1.0, 1.0], [-10.0, 12.0])
```

If there is only one partition (like above) for each variable, then the function formulates the McCormick relaxation. The package supports more than one partition only one of the variables `x` or `y`. Support for multi-variable partitioning will be added in the future versions. The following example illustrates the API for more than one partitions on the `x2` variable:

```julia 
using PolyhderalRelaxations, JuMP
m = Model() 
@variable(m, -1.0 <= x1 <= 1.0)
@variable(m, -1.0 <= x2 <= 1.0)
@variable(m, z)
construct_bilinear_relaxation!(m, x1, x2, z, [-1.0, 1.0], [-1.0, -0.25, 0.25, 1.0])
```

For examples, the reader is referred to the unit tests in the `test/` folder of the GitHub repository. 

## API for the MILP/LP Relaxation for Nonlinear Univariate function
The API documentation for both the MILP and the LP relaxations for a given nonlinear univariate function is as follows:

```@docs 
construct_univariate_relaxation!
```

The derivate of the function is an optional keyword argument. If the derivate is not specified, the package invokes `ForwardDiff.jl` to compute the derivative.

## API for the MILP/LP Relaxation for Bilinear Term
The API for both the MILP and the LP relaxations for a bilinear term is as follows:

```@docs 
construct_bilinear_relaxation!
```

## Additional Algorithms for Relaxations of Univariate Functions
Apart from providing MILP and LP relaxations of graphs of univariate function ``y=f(x)`` for a given partition, the package also implements some basic partitioning algorithms to refine the partition and provide tighter relaxations. Given a partition, PolyhedralRelaxations.jl can refine the partition using an interval-bisection algorithm detailed in the following reference: 

* K. Sundar, S. Sanjeevi, and H, Nagarajan (2020). Sequence of Polyhedral Relaxations for Nonlinear Univariate Functions. ([arxiv link](https://arxiv.org/abs/2005.13445))

The partition refinement scheme is equipped with multiple stopping criteria that can be toggled by setting non-default values to the following keyword arguments in the standard API:

1. `error_tolerance`: this parameter measures the maximum vertical distance between the under- and over-estimators of the relaxation within the domain of the univariate function. By default, this value is set to `NaN64`. Setting a non-zero positive value for this parameter would result partition and a corresponding relaxation such that the maximum vertical distance between the under- and over-estimators of the relaxation is at most the value prescribed by `error_tolerance`. 

2. `num_additional_partitions`: this parameter as the name suggests will refine the partition using the interval bisection algorithm till the total number of additional partitions reaches this number. It's default value is set to 0. 

The refinement algorithm will stop if either of the above two stopping criteria is satisfied. 
The following example illustrates the API for refining the partition using interval-bisection for a univariate function. 

```julia 
using PolyhderalRelaxations, 
partition = [1.0, 2.0]
refine_partition!(x -> x^2, partition, num_additional_partitions=2)
```

The new partition would now be `[0.0, 0.5, 0.75, 1.0]`.

## API for the partition refinement
The API for partition refinement using interval bisection is as follows:

```@docs 
refine_partition!
```

## Tolerance parameters
The main functions to produce the relaxations are also equipped with the following two tolerance parameters:

1. `length_tolerance`: this parameter is used to control the refinement algorithm. A partition whose length is less than the `length_tolerance` is not partitioned further. The default value of this parameter is ``\epsilon = 1 \times 10^{-6}``.

2. `derivative_tolerance`: when the refinement algorithm is used with the `base_partition` not containing the inflection points of the function in its domain, the resulting relaxations will be erroneous i.e., there is no guarantee that the relaxation obtained is even valid. One necessary condition to detect this issue is by ensuring that the derivatives at successive partition points are not equal. This condition is checked up to the tolerance specified by this parameter. The default value of this parameter is ``\epsilon = 1 \times 10^{-6}``.






