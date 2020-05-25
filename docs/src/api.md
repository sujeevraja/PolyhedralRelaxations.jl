API
=================


## Basic Usage
Here is an example illustrating the basic API 

```julia 
using PolyhderalRelaxations
milp_relaxation, function_data = construct_milp_relaxation(x -> x^3, [-1.0, 0.0, 1.0])
lp_relaxation, function_data = construct_lp_relaxation(x -> x^2, [-1.0, 0.0, 1.0])
```

For more examples and details on the helper functions, the reader is referred to the unit tests in the `test/` folder of the GitHub repository.


## Additional Algorithms
Apart from providing MILP and LP relaxations of graphs of univariate function ``y=f(x)`` for a given partition, the package also implements some basic partitioning algorithms to refine the partition and provide tighter relaxations. Given a partition, PolyhedralRelaxations.jl can refine the partition using an interval-bisection algorithm detailed in the following reference: 

* K. Sundar, S. Sanjeevi, and H, Nagarajan (forthcoming). Sequence of Polyhedral Relaxations for Nonlinear Univariate Functions.

The partition refinement scheme is equipped with multiple stopping criteria that can be toggled by setting non-default values to the following keyword arguments in the standard API:

1. `error_tolerance`: this parameter measures the maximum vertical distance between the under- and over-estimators of the relaxation within the domain of the univariate function. By default, this value is set to `NaN64`. Setting a non-zero positive value for this parameter would result partition and a corresponding relaxation such that the maximum vertical distance between the under- and over-estimators of the relaxation is at most the value prescribed by `error_tolerance`. 

2. `num_additional_partitions`: this parameter as the name suggests will refine the partition using the interval bisection algorithm till the total number of additional partitions reaches this number. It's default value is set to 0. 

The refinement algorithm will stop if either of the above two stopping criteria is satisfied. 

## Tolerance parameters
The main functions to produce the relaxations are also equipped with the following two tolerance parameters:

1. `length_tolerance`: this parameter is used to control the refinement algorithm. A partition whose length is less than the `length_tolerance` is not partitioned further. The default value of this parameter is ``\epsilon = 0.001``.

2. `derivative_tolerance`: when the refinement algorithm is used with the `base_partition` not containing the inflection points of the function in its domain, the resulting relaxations will be erroneous i.e., there is no guarantee that the relaxation obtained is even valid. One necessary condition to detect this issue is by ensuring that the derivatives at successive partition points are not equal. This condition is checked up to the tolerance specified by this parameter. The default value of this parameter is ``\epsilon = 0.001``.

## API for the MILP Relaxation
The two main functions to obtain the MILP relaxation of a given univariate function are 

```@docs 
construct_milp_relaxation
```

The function without the derivative invokes the function with the derivative after computing it using the Julia package `ForwardDiff.jl`

## API for the LP Relaxation
The main functions to obtain the LP relaxation of a given univariate function are 

```@docs 
construct_lp_relaxation
```
The function without the derivative invokes the function with the derivative after computing it using the Julia package `ForwardDiff.jl`



