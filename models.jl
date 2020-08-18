using Distributions

abstract type Models end


function find_x(model::Models, cones::AbstractVector{Cone}, pars::Pars; kwargs...)

    scores   = scores(model, cones, pars; kwargs...)
    ~, i_max = findmax(scores)

    return deepcopy(cones[i_max].center)
end


mean_folded_normal(μ::Real, σ_sq::Real) = sqrt(2*σ_sq/pi)*exp(-μ^2/(2*σ_sq)) + μ*(1-2*cdf(Normal(0,1), -μ/(sqrt(σ_sq))))


function_variance(cones, pars) = mean_folded_normal.(pars.s.*width(cones), pars.h_s.^2 .+ pars.s.^2 .* pars.w_s)
height_variance(cones, pars)   = mean_folded_normal.(0, pars.h_s.^2)


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


scores(model::Robust1, cones::AbstractVector{Cone}, pars::Pars; δ::Real, kwargs...) = (height(cones).-δ)./abs.(pars.s .* width(cones))


################
# Robust 2
################


struct Robust2 <: Models end


scores(model::Robust2, cones::AbstractVector{Cone}, pars::Pars; δ::Real, kwargs...) = (height(cones).-δ)./height_variance(cones, pars)


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


scores(model::Yazdani1, cones::AbstractVector{Cone}, pars::Pars; kwargs...) = height(cones) .- height_variance(cones, pars)


################
# Yazdani 2
################


struct Yazdani2 <: Models end


scores(model::Yazdani2, cones::AbstractVector{Cone}, pars::Pars; kwargs...) = -pars.s


################
# Yazdani 3
################


struct Yazdani3 <: Models end


scores(model::Yazdani3, cones::AbstractVector{Cone}, pars::Pars; kwargs...) = .-height_variance(cones, pars)


################
# Yazdani 4
################


struct Yazdani4 <: Models end


function scores(model::Yazdani4, cones::AbstractVector{Cone}, pars::Pars; kwargs...)

    hv = height_variance(cones, pars)
    return -(hv/maximum(hv) .+ pars.s/maximum(pars.s))
end

################
# Yazdani 5
################


struct Yazdani5 <: Models end


function scores(model::Yazdani5, cones::AbstractVector{Cone}, pars::Pars; δ::Real, kwargs...)

    s       = zeros(pars.m)
    ii      = height(cones) .>= δ
    hv      = height_variance(cones, pars)
    s[ii]   = -(hv[ii]./maximum(hv) .+ pars.s[ii]./maximum(pars.s))
    s[.!ii] = hv[.!ii] .- maximum(height(cones)) .- 2.
    return s
end

################
# Utilities
################


maxk_i(a, k) = partialsortperm(a, 1:k, rev=true)
