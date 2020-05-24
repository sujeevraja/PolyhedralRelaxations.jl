Quick Start Guide
=================

## Geometry of relaxations
This quick start guide will introduce the main concepts of PolyhedralRelaxations. The package constructs MILP and LP relaxations for ``y = f(x)``, given the function and discretization points on the domain of ``f(x)``. For instance, given ``y=x^3`` and a partition of the domain ``[-1.5, 0, 2]``, it constructs an MILP relaxation of ``y = x^3`` as the disjunction of the triangles shown in the figure below:

![example](assets/example.svg)

Notice from the above example that, for the above disjunction of triangles to be a valid relaxation for ``y = f(x)``, the partition points must necessarily include the inflection points of the univariate function. Also, in the above example, the domain of the function ``y = x^3`` is given by ``[-1.5, 2]``. The following table provides a list of inflection points for various well-used univariate functions. 

| Function | Inflection points |
| --- | --- |
| ``x^n`` (``n`` is odd)  | ``\{0\}`` |
| ``x^n`` (``n`` is even) | ``\phi`` |
| ``\sin x`` | ``n\pi`` ``(n \in \mathbb Z)`` | 
| ``\cos x`` | ``n\frac{\pi}{2}`` ``(n`` is odd ``)``| 
| ``x\|x\|`` | ``\{0\}``| 

If the domain of the provided function contains an inflection point, then the requirement is that the inflection points should be provided as a discretization point; failing which the MILP relaxation need not be valid for the given function.

As far as the LP relaxation is concerned, the package will return the convex hull of the triangles as a linear (sparse) constraint matrix and an associated RHS vector. 

