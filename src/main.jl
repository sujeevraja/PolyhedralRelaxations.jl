using Memento
using SparseArrays

# logger_config!("debug")

"""
Evaluate the value of the given univariate function `f` at each point `x` in
the given set of partition points and return a list of vertices ((x,f(x))) as
Vertex objects.
"""
function collect_secant_vertices(
    f::Function, partition_points::Array{Real,1})::Array{Vertex,1}
    secant_vertices = Vertex[]
    num_points = length(partition_points)
    for x in partition_points
        push!(secant_vertices, Vertex(x, f(x)))
        v = secant_vertices[end]
        info(_LOGGER, "x: $x v_x: $(v.x) v_y: $(v.y) ")
    end
    return secant_vertices
end

"""
Use the given derivative function `f_dash` to generate and return tangent
vertices as a list of Vertex objects.
"""
function collect_tangent_vertices(
    f_dash::Function,
    secant_vertices::Array{Vertex,1})::Array{Vertex,1}
    # tangent_vertices[i] is the vertex formed by the intersection of tangents
    # of the given curve f at the partition interval formed by thex
    # x-coordinates of secant_vertices[i] and secant_vertices[i+1].
    # In terms of the notation in the paper, if secant_vertices[i] is the
    # vertex v_i and secant_vertices[i+1] is the vertex v_{i+1},
    # tangent_vertices[i] is the vertex v_{i,i+1}.
    num_points = length(secant_vertices)
    tangent_vertices = Vertex[]

    for i in 1:num_points-1
        v0 = secant_vertices[i]
        v1 = secant_vertices[i+1]

        # Find derivative values at secant points. By definition, the partition
        # cannot have
        d0 = f_dash(v0.x)
        d1 = f_dash(v1.x)
        @assert !isapprox(d0, d1, atol=1e-5)

        # Compute x-coordinate of tangent vertex. This is the intersection of
        # the tangents to the curve at v0.x and v1.x.
        x_new = (v1.y - v0.y + (d0*v0.x) - (d1*v1.x)) / (d0 - d1)
        @assert x_new >= v0.x
        @assert x_new <= v1.x

        y_new = v0.y + (d0*(x_new - v0.x))
        push!(tangent_vertices, Vertex(x_new, y_new))
        info(_LOGGER,  "x0: $(v0.x) x1: $(v1.x) x_new: $x_new y_new: $y_new")
    end
    return tangent_vertices
end

function build_model(
        uf::UnivariateFunction,
        secant_vertices::Vector{Vertex},
        tangent_vertices::Vector{Vertex})::Model
    info(_LOGGER, "starting to build model...")

    # Indices to recover variable values from model. Indices of delta_1^i,
    # delta_2^i and z_i start from 1.
    x_index = 1
    y_index = 2
    num_points = length(secant_vertices)
    info(_LOGGER, "number of partition points: $num_points")

    start = 3
    delta_1_indices = collect(start:(num_points+start))
    info(_LOGGER, "delta_1_indices: $start to $(delta_1_indices[end])")

    start = delta_1_indices[end]+1
    delta_2_indices = collect(start:(num_points+start))
    info(_LOGGER, "delta_2_indices: $start to $(delta_2_indices[end])")

    start = delta_2_indices[end]+1
    z_indices = collect(start:(num_points+start))
    info(_LOGGER, "z_indices: $start to $(z_indices[end])")

    # Constraint coefficient data
    row_indices = Int64[]
    col_indices = Int64[]
    coefs = Real[]

    # RHS data
    rhs_row_indices = Int64[]
    rhs_values = Real[]

    # constraint rows are ordered as:
    # x constraint
    # y constraint
    # \delta_1^1 + \delta_2^1 <= 1
    # \delta_1^i + \delta_2^i <= z_{i-1}
    # z_{i-1} <= \delta_2^{i-1}

    # ------ x constraint data ------
    row = 1
    push!(row_indices, row)
    push!(col_indices, x_index)
    # -------------------------------

    # Store constraint data into a Model object and return it.
    A = sparse(row_indices, col_indices, coefs)
    b = sparsevec(rhs_row_indices, rhs_values)
    info(_LOGGER, "completed building model.")

    return Model(A, b, x_index, y_index,
        delta_1_indices, delta_2_indices, z_indices)
end

function main()
    info(_LOGGER, "starting model generation...")

    lb = -2.0
    ub = 2.0
    uf = UnivariateFunction(
        x->x^3,  # f
        x->3 * (x^2),  # f'
        domain_lb = lb,
        domain_ub = ub,
        inflection_points = Array{Real}(collect(lb:0.25:ub)))

    sec_vertices = collect_secant_vertices(uf.f, uf.inflection_points)
    info(_LOGGER, "collected $(length(sec_vertices)) secant vertices.")

    tan_vertices = collect_tangent_vertices(uf.f_dash, sec_vertices)
    info(_LOGGER, "collected $(length(tan_vertices)) tangent vertices.")

    # model = build_model(uf, secant_vertices, tangent_vertices)
    info(_LOGGER, "completed model generation.")
end
