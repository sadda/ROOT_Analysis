abstract type Models end


################
# TMO
################


struct TMO <: Models end


function find_x(model::TMO, cones::AbstractVector{Cone}, pars::Pars; kwargs...)

    ~, i_max = findmax(height(cones))
    x = cones[i_max].center

    return deepcopy(x)
end


################
# Robust 1
################


struct Robust1 <: Models end


function find_x(model::Robust1, cones::AbstractVector{Cone}, pars::Pars; δ::Real, kwargs...)

    ~, i_max = findmax((height(cones).-δ)./abs.(width(cones)))
    x = cones[i_max].center

    return deepcopy(x)
end


################
# RoA 1
################


struct RoA1 <: Models

    n_eval::Int  # First round of evaluation
    n_sel::Int   # Number of points for improvement in the second round

    function RoA1(n::Int)

        n_eval = round(n*0.8)
        n_sel  = 5
        return new(n_eval, n_sel)
    end
end


update_x(model::RoA1, cones::AbstractVector{Cone}, x) = cones[cones_max_i(x, cones)].center


function find_x(model::RoA1, cones::AbstractVector{Cone}, pars::Pars; kwargs...)

    x_try = random_points(model.n_eval, pars)
    f_try = cones_value(x_try, cones)
    i_try = maxk_i(f_try, model.n_sel)

    x_red = x_try[:,i_try]

    x_upd = hcat([update_x(model, cones, x) for x in eachcol(x_red)]...)
    f_upd = cones_value(x_upd, cones)

    ~, i_max = findmax(f_upd)
    x        = x_upd[:,i_max]

    return deepcopy(x)
end


################
# Yazdani 1
################


struct Yazdani1 <: Models end


function find_x(model::Yazdani1, cones::AbstractVector{Cone}, pars::Pars; kwargs...)

    ~, i_max = findmax(height(cones) .- sqrt(2/pi)*pars.h_s)
    x        = cones[i_max].center

    return deepcopy(x)
end


################
# Yazdani 2
################


struct Yazdani2 <: Models end


function find_x(model::Yazdani2, cones::AbstractVector{Cone}, pars::Pars; kwargs...)

    ~, i_min = findmin(pars.s)
    x        = cones[i_min].center

    return deepcopy(x)
end


################
# Yazdani 3
################


struct Yazdani3 <: Models end


function find_x(model::Yazdani3, cones::AbstractVector{Cone}, pars::Pars; kwargs...)

    ~, i_min = findmin(pars.h_s)
    x        = cones[i_min].center

    return deepcopy(x)
end


################
# Yazdani 4
################


struct Yazdani4 <: Models end


function find_x(model::Yazdani4, cones::AbstractVector{Cone}, pars::Pars; kwargs...)

    ~, i_min = findmin(pars.h_s/maximum(pars.h_s) .+ pars.s/maximum(pars.s))
    x        = cones[i_min].center

    return deepcopy(x)
end


################
# Yazdani 5
################


struct Yazdani5 <: Models end


function find_x(model::Yazdani5, cones::AbstractVector{Cone}, pars::Pars; δ::Real, kwargs...)

    ii       = findall(height(cones) .>= δ)
    if length(ii) > 0
        ~, i_min = findmin(pars.h_s[ii]./maximum(pars.h_s[ii]) .+ pars.s[ii]./maximum(pars.s[ii]))
        x        = cones[ii[i_min]].center
    else # All cones are too small
        ~, i_max = findmax(height(cones) .- sqrt(2/pi)*pars.h_s)
        x        = cones[i_max].center
    end

    return deepcopy(x)
end


################
# Utilities
################


maxk_i(a, k) = partialsortperm(a, 1:k, rev=true)
