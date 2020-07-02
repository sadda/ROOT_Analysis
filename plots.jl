using Plots


function plot_contours(cones::AbstractVector{Cone}, pars::Pars; δ_x=0.25, min_value=-Inf, plot_boundary=false)

    pars.d != 2 && error("Only two dimensional plotting supported.")

    x = pars.x_min[1]:δ_x:pars.x_max[1]
    y = pars.x_min[2]:δ_x:pars.x_max[2]
    n_x = length(x)
    n_y = length(y)

    f(x, y) = max(cones_value([x, y], cones), min_value)

    X = repeat(reshape(x, 1, :), n_y, 1)
    Y = repeat(y, 1, n_x)
    Z = map(f, X, Y)

    p = contour(x, y, Z, levels = 20, c=:jet)

    if plot_boundary

        II = [cones_max_i([X[i,j], Y[i,j]], cones) for i=1:n_x, j=1:n_x]
        JJ = [II[i,j]==II[i-1,j] && II[i,j]==II[i+1,j] && II[i,j]==II[i,j-1] && II[i,j]==II[i,j+1] for i=2:n_x-1, j=2:n_x-1]

        KK = falses(size(II))
        for i=1:n_x, j=1:n_x
            cond1 = i==1 || i==n_x || j==1 || j==n_x
            if cond1 || II[i,j]!=II[i-1,j] || II[i,j]!=II[i+1,j] || II[i,j]!=II[i,j-1] || II[i,j]!=II[i,j+1]
                KK[i,j] = true
            end
        end
        scatter!(X[KK], Y[KK], markersize=3, markercolor=:black, legend=:none)
    end

    return p
end


function add_peaks(p::AbstractPlot, cones::AbstractVector{Cone}, pars::Pars; kwargs...)

    centers = center(cones)
    scatter!(p, centers[1,:], centers[2,:]; markershape=:cross, markercolor=:black, kwargs...)
end
