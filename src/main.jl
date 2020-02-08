using Memento
using SparseArrays

# logger_config!("debug")

"""
    collect_secant_vertices(f, partition_points)

Evaluate the value of the univariate function `f` at each point x in
`partition_points` and return a list of vertices of the form (x, f(x)) as
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
    collect_tangent_vertices(f_dash, secant_vertices)

Generate tangent vertices for each interval formed by x-coordinates of
adjacent vertices in `secant_vertices`. Then, return the generated tangent
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

"""
    build_model(uf, secant_vertices, tangent_vertices)

Collect constraint data of the MILP formulation of the polyhedral relaxation
in a Model object and return it.
"""
function build_model(
        uf::UnivariateFunction,
        secant_vertices::Vector{Vertex},
        tangent_vertices::Vector{Vertex})::Model
    info(_LOGGER, "starting to build model...")

    # Indices to recover variable values from model. Indices of delta_1^i,
    # delta_2^i and z_i start from 1.
    num_points = length(secant_vertices)
    index_data = IndexData(num_points)
    info(_LOGGER, "number of partition points: $num_points")

    i_start = index_data.delta_1_indices[1]
    i_end = index_data.delta_1_indices[end]
    info(_LOGGER, "delta_1_indices: $i_start to $i_end")

    i_start = index_data.delta_2_indices[1]
    i_end = index_data.delta_2_indices[end]
    info(_LOGGER, "delta_2_indices: $i_start to $i_end")

    i_start = index_data.z_indices[1]
    i_end = index_data.z_indices[end]
    info(_LOGGER, "z_indices: $i_start to $i_end")

    constraint_data = ConstraintData()
    add_x_constraint(constraint_data, index_data, secant_vertices,
        tangent_vertices)
    add_y_constraint(constraint_data, index_data, secant_vertices,
        tangent_vertices)
    add_first_delta_constraint(constraint_data, index_data)

    # constraint rows are ordered as:
    # \delta_1^1 + \delta_2^1 <= 1
    # \delta_1^i + \delta_2^i <= z_{i-1}
    # z_{i-1} <= \delta_2^{i-1}

    # Ensure that constraint symbols are valid.
    for s in constraint_data.constraint_senses:
        @assert s in possible_senses
    end

    # Store constraint data into a Model object and return it.
    A = sparse(constraint_data.constraint_row_indices,
        constraint_data.constraint_column_indices,
        constraint_data.constraint_coefficients)
    b = sparsevec(constraint_data.rhs_row_indices, constraint_data.rhs_values)
    info(_LOGGER, "completed building model.")

    return Model(A, b,
        index_data.x_index,
        index_data.y_index,
        index_data.delta_1_indices,
        index_data.delta_2_indices,
        index_data.z_indices,
        constraint_data.constraint_senses)
end

function add_x_constraint(
        constraint_data::ConstraintData,
        index_data::IndexData,
        secant_vertices::Vector{Vertex},
        tangent_vertices::Vector{Vertex})
    row = constraint_data.num_constraints+1

    # Add x variable to constraint.
    add_coef(constraint_data, row, index_data.x_index, 1.0)

    num_vars = length(secant_vertices) - 1
    for i in 1:num_vars
        # Add delta_1 variable to constraint.
        column = index_data.delta_1_indices[i]
        value = secant_vertices[i].x - tangent_vertices[i].x
        add_coef(constraint_data, row, column, value)

        # Add delta_2 variable to constraint.
        column = index_data.delta_2_indices[i]
        value = secant_vertices[i].x - secant_vertices[i+1].x
        add_coef(constraint_data, row, column, value)
    end

    # Add sense
    push!(constraint_data.constraint_senses, :eq)

    # Add right hand side.
    add_rhs(constraint_data, row, secant_vertices[1].x)

    constraint_data.num_constraints += 1
    info(_LOGGER, "built x coordinate constraint.")
end


function add_y_constraint(
        constraint_data::ConstraintData,
        index_data::IndexData,
        secant_vertices::Vector{Vertex},
        tangent_vertices::Vector{Vertex})
    row = constraint_data.num_constraints+1

    # Add y variable to constraint.
    add_coef(constraint_data, row, index_data.y_index, 1.0)

    num_vars = length(secant_vertices) - 1
    for i in 1:num_vars
        # Add delta_1 variable to constraint.
        column = index_data.delta_1_indices[i]
        value = secant_vertices[i].y - tangent_vertices[i].y
        add_coef(constraint_data, row, column, value)

        # Add delta_2 variable to constraint.
        column = index_data.delta_2_indices[i]
        value = secant_vertices[i].y - secant_vertices[i+1].y
        add_coef(constraint_data, row, column, value)
    end

    # Add sense
    push!(constraint_data.constraint_senses, :eq)

    # Add right hand side.
    add_rhs(constraint_data, row, secant_vertices[1].y)

    constraint_data.num_constraints += 1
    info(_LOGGER, "built y coordinate constraint.")
end

function add_first_delta_constraint(
        constraint_data::ConstraintData,
        index_data::IndexData)
    row = constraint_data.num_constraints + 1
    add_coef(constraint_data, row, index_data.delta_1_indices[1], 1.0)
    add_coef(constraint_data, row, index_data.delta_2_indices[1], 1.0)
    push!(constraint_data.constraint_senses, :leq)
    add_rhs(constraint_data, row, 1.0)
    constraint_data.num_constraints += 1
    info(_LOGGER, "built delta_1^1 + delta_2^1 <= 1 constraint.")
end

function add_coef(
        constraint_data::ConstraintData, row::Int64, col::Int64, value::Real)
    push!(constraint_data.constraint_row_indices, row)
    push!(constraint_data.constraint_column_indices, col)
    push!(constraint_data.constraint_coefficients, value)
end

function add_rhs(constraint_data::ConstraintData, row::Int64, value::Real)
    push!(constraint_data.rhs_row_indices, row)
    push!(constraint_data.rhs_values, value)
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

    model = build_model(uf, sec_vertices, tan_vertices)
    info(_LOGGER, "completed model generation.")
end
