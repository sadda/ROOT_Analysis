abstract type Models end


function find_x(model::Models, cones::AbstractVector{Cone}, pars::Pars; kwargs...)

    scores   = scores(model, cones, pars; kwargs...)
    ~, i_max = findmax(scores)

    return deepcopy(cones[i_max].center)
end


################
# TMO
################


struct TMO <: Models end


scores(model::TMO, cones::AbstractVector{Cone}, pars::Pars; kwargs...) = height(cones)


################
# Random
################


struct Random <: Models end


scores(model::Random, cones::AbstractVector{Cone}, pars::Pars; kwargs...) = rand(pars.m)


################
# Random above
################


struct RandomAbove <: Models end


function scores(model::RandomAbove, cones::AbstractVector{Cone}, pars::Pars; δ::Real, kwargs...)

    s       = zeros(pars.m)
    ii      = height(cones) .>= δ
    s[ii]   = rand(sum(ii))
    s[.!ii] = -rand(sum(.!ii))
    return s
end


################
# Robust 1
################


struct Robust1 <: Models end


scores(model::Robust1, cones::AbstractVector{Cone}, pars::Pars; δ::Real, kwargs...) = (height(cones).-δ)./abs.(width(cones))


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


function scores(model::RoA1, cones::AbstractVector{Cone}, pars::Pars; n=1000, kwargs...)

    x      = random_points(n, pars)
    i_cone = cones_max_i(x, cones)
    s      = [sum(i_cone .== i) for i in 1:pars.m]
    return s./sum(s)
end


################
# Yazdani 1
################


struct Yazdani1 <: Models end


scores(model::Yazdani1, cones::AbstractVector{Cone}, pars::Pars; kwargs...) = height(cones) .- sqrt(2/pi)*pars.h_s


################
# Yazdani 2
################


struct Yazdani2 <: Models end


scores(model::Yazdani2, cones::AbstractVector{Cone}, pars::Pars; kwargs...) = -pars.s


################
# Yazdani 3
################


struct Yazdani3 <: Models end


scores(model::Yazdani3, cones::AbstractVector{Cone}, pars::Pars; kwargs...) = -pars.h_s


################
# Yazdani 4
################


struct Yazdani4 <: Models end


scores(model::Yazdani4, cones::AbstractVector{Cone}, pars::Pars; kwargs...) = -(pars.h_s/maximum(pars.h_s) .+ pars.s/maximum(pars.s))


################
# Yazdani 5
################


struct Yazdani5 <: Models end


function scores(model::Yazdani5, cones::AbstractVector{Cone}, pars::Pars; δ::Real, kwargs...)

    s       = zeros(pars.m)
    ii      = height(cones) .>= δ
    s[ii]   = -(pars.h_s[ii]./maximum(pars.h_s[ii]) .+ pars.s[ii]./maximum(pars.s[ii]))
    s[.!ii] = (height(cones)[.!ii] .- sqrt(2/pi)*pars.h_s[.!ii]) .- maximum(height(cones)) .- 2.
    return s
end

################
# Utilities
################


maxk_i(a, k) = partialsortperm(a, 1:k, rev=true)
