PolyhedralRelaxations.jl Change Log
=========================

### v0.4.0
- Add dependabot support
- Fix bug due to lack of ``JuMP.is_fixed`` support
- Add multilinear LP and MILP (SOS-2) relaxations 
- Add multilinear linking constraints given multiple terms 
- Add partition refinement schemes (:bisect_all, :bisect, :at_point, :non_uniform)
- Add formulation reuse for univariate and bilinear relaxations

### v0.3.5 
- Add compat entries for test dependencies
- Move Cbc test dependency to HiGHS
- Allow for JuMP v1.X 
- Update versioning for LoggingExtras.jl

### v0.3.4
- Allow for JuMP 1.0
- Move logging to Logging.jl and LoggingExtras.jl from Memento
- Fix tests
- Drop Memento from dependencies
- Add linting and formatting workflows

### v0.3.3
- Allow for JuMP 0.23

### v0.3.2
- Allow for JuMP 0.22

### v0.3.1
- Documentation updates
- Adds base name prefixes as optional argument
- Memento version update

### v0.3.0
- Add piecewise bilinear relaxations with partition on one variable
- API support for the above

### v0.2.2
- Bug fix in documentation generator

### v0.2.1
- Example bug fix in LP relaxation

### v0.2.0 
- API refactor 
- Move CI to GitHub Actions 
- Make JuMP a dependency and populate relaxations in JuMP model directly 
- Limit exports
- Clean up documentation for new API
- Refactor examples

### v0.1.2
- patch release
- Update DataStructures version 

### v0.1.1
- patch release
- Update Memento version
- Clean up citation in readme and documentation

### v0.1.0
- Initial implementation 
