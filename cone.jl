using Distributions
using LinearAlgebra


mutable struct Cone

    d::Int
    center::AbstractVector
    height::Real
    width::Real

    function Cone(center::AbstractVector{<:Float64}, height::Float64, width::Float64)

        d = length(center)
        return new(d, center, height, width)
    end
end


center(cone::Cone) = cone.center
height(cone::Cone) = cone.height
width(cone::Cone)  = cone.width


center(cones::AbstractVector{Cone}) = hcat(center.(cones)...)
height(cones::AbstractVector{Cone}) = height.(cones)
width(cones::AbstractVector{Cone})  = width.(cones)


function init_cones(pars::Pars)

    centers = pars.x_min .+ (pars.x_max-pars.x_min).*rand(pars.d,pars.m)
    heights = pars.h_init_min .+ (pars.h_init_max-pars.h_init_min).*rand(pars.m)
    widths  = pars.w_init_min .+ (pars.w_init_max-pars.w_init_min).*rand(pars.m)

    a = [Cone(centers[:,k], heights[k], widths[k]) for k in 1:m]
end


function update_reflect!(x::V, step::V, x_min::V, x_max::V) where {T<:Real, V<:AbstractVector{T}}

    x .+= step

    ii    = x .< x_min
    x[ii] .+= 2*(x_min[ii] - x[ii])

    ii    = x .> x_max
    x[ii] .-= 2*(x[ii] - x_max[ii])
end


function update_reflect(x::T, step::T, x_min::T, x_max::T) where {T<:Real}

    x += step

    if x < x_min
        x += 2*(x_min - x)
    elseif x > x_max
        x -= 2*(x - x_max)
    end

    return x
end


function update_cone!(cone::Cone, pars::Pars, s::Real, h_s::Real, w_s::Real)

    step_x = s*pars.step_x_generator()
    step_h = h_s*pars.step_h_generator()
    step_w = w_s*pars.step_w_generator()

    update_reflect!(cone.center, step_x, pars.x_min, pars.x_max)
    cone.height = update_reflect(cone.height, step_h, pars.h_min, pars.h_max)
    cone.width  = update_reflect(cone.width, step_w, pars.w_min, pars.w_max)
end


function update_cones!(cones::AbstractVector{Cone}, pars::Pars)

    return update_cone!.(cones, [pars], pars.s, pars.h_s, pars.w_s)
end


function cone_value(x::AbstractVector, cone::Cone)

    return cone.height - cone.width*norm(x-cone.center)
end


function cones_value(x::AbstractVector, cones::AbstractVector{Cone})

    return maximum(cone_value.([x], cones))
end


function cones_value(x_mat::AbstractMatrix, cones::AbstractVector{Cone})

    return [cones_value(x, cones) for x in eachcol(x_mat)]
end


function cones_max_i(x::AbstractVector, cones::AbstractVector{Cone})

    ~, i_max = findmax(cone_value.([x], cones))
    return i_max
end


function cones_max_i(x_mat::AbstractMatrix, cones::AbstractVector{Cone})

    return [cones_max_i(x, cones) for x in eachcol(x_mat)]
end
