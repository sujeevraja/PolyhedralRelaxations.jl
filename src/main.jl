using Memento

# logger_config!("debug")

function collect_vertices(uf::UnivariateFunction)
    info(_LOGGER, "collecting secant vertices...")
    secant_vertices = Vertex[]
    num_points = length(uf.inflection_points)
    for x in uf.inflection_points
        push!(secant_vertices, Vertex(x, uf.f(x)))
        v = secant_vertices[end]
        info(_LOGGER, "x: $x v_x: $(v.x) v_y: $(v.y) ")
    end
    info(_LOGGER, "collected secant_vertices.")

    info(_LOGGER, "collecting tangent vertices...")
    # tangent_vertices[i] is the vertex formed by the intersection of tangents
    # of the given curve f at the partition interval formed by thex
    # x-coordinates of secant_vertices[i] and secant_vertices[i+1].
    # In terms of the notation in the paper, if secant_vertices[i] is the
    # vertex v_i and secant_vertices[i+1] is the vertex v_{i+1},
    # tangent_vertices[i] is the vertex v_{i,i+1}.
    tangent_vertices = Vertex[]
    for i in 1:num_points-1
        v0 = secant_vertices[i]
        v1 = secant_vertices[i+1]

        # Find derivative values at secant points. By definition, the partition
        # cannot have
        d0 = uf.f_dash(v0.x)
        d1 = uf.f_dash(v1.x)
        @assert !isapprox(d0, d1, atol=1e-5)

        # Compute x-coordinate of tangent vertex. This is the intersection of
        # the tangents to the curve at v0.x and v1.x.
        x_new = (v1.y - v0.y + (d0*v0.x) - (d1*v1.x)) / (d0 - d1)
        @assert x_new >= v0.x
        @assert x_new <= v1.x

        y_new = v0.y + (d0*(x_new - v0.x))
        curve_new = uf.f(x_new)
        push!(tangent_vertices, Vertex(x_new, y_new))
        info(_LOGGER, "x0: $(v0.x) x1: $(v1.x) x_new: $x_new y_new: $y_new, curve_new: $curve_new")
    end
    info(_LOGGER, "collected tangent vertices.")

    return Pair(secant_vertices, tangent_vertices)
end

function main()
    info(_LOGGER, "starting model generation...")

    uf = UnivariateFunction(
        x->x^3,  # f
        x->3 * (x^2),  # f'
        domain_lb = -2.0,
        domain_ub = 2.0,
        inflection_points = Array{Real}(collect(-2:0.25:2)))

    (secant_vertices, tangent_vertices) = collect_vertices(uf)
    info(_LOGGER, "completed model generation.")
end
