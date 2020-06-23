struct Pars{T<:Real}
    d::Int
    m::Int
    x_min::AbstractVector{T}
    x_max::AbstractVector{T}
    h_min::T
    h_max::T
    w_min::T
    w_max::T
    h_init_min::T
    h_init_max::T
    w_init_min::T
    w_init_max::T
    s::AbstractVector{T}
    h_s::AbstractVector{T}
    w_s::AbstractVector{T}

    function Pars(d::Int, m::Int, x_min, x_max, h_min, h_max, w_min, w_max, h_init, w_init, s_min, s_max, h_s_min, h_s_max, w_s_min, w_s_max)

        x_min = repeat_n(x_min, d)
        x_max = repeat_n(x_max, d)

        h_init_min, h_init_max = modify_init(h_init, h_min, h_max)
        w_init_min, w_init_max = modify_init(w_init, w_min, w_max)

        s_min = repeat_n(s_min, m)
        s_max = repeat_n(s_max, m)
        s     = s_min .+ (s_max-s_min).*rand(m)

        h_s_min = repeat_n(h_s_min, m)
        h_s_max = repeat_n(h_s_max, m)
        h_s     = h_s_min .+ (h_s_max-h_s_min).*rand(m)

        w_s_min = repeat_n(w_s_min, m)
        w_s_max = repeat_n(w_s_max, m)
        w_s     = w_s_min .+ (w_s_max-w_s_min).*rand(m)

        @assert length(x_min) == d
        @assert length(x_max) == d
        @assert length(s)     == m
        @assert length(h_s)   == m
        @assert length(w_s)   == m

        return new{typeof(h_min)}(d, m, x_min, x_max, h_min, h_max, w_min, w_max, h_init_min, h_init_max, w_init_min, w_init_max, s, h_s, w_s)
    end
end


function modify_init(h_init, h_min, h_max)
    if isa(h_init, Number)
        h_init_min = h_init
        h_init_max = h_init
    elseif h_init == "uniform"
        h_init_min = h_min
        h_init_max = h_max
    else
        error("Init value not defined well")
    end
    return h_init_min, h_init_max
end


repeat_n(x::Real, n) = x*ones(n)
repeat_n(x::AbstractVector, n) = x
