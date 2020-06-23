abstract type Models end


################
# TMO
################


struct TMO <: Models end


function find_x(model::TMO, cones::AbstractVector{Cone}, pars::Pars)

    ~, i_max = findmax(height(cones))
    x = cones[i_max].center

    return deepcopy(x)
end


################
# Robust 1
################


struct Rob1 <: Models

    n_eval::Int  # First round of evaluation
    n_sel::Int   # Number of points for improvement in the second round

    function Rob1(n::Int)

        n_eval = round(n*0.8)
        n_sel  = 5
        return new(n_eval, n_sel)
    end
end


update_x(model::Rob1, cones::AbstractVector{Cone}, x) = cones[cones_max_i(x, cones)].center


function find_x(model::Rob1, cones::AbstractVector{Cone}, pars::Pars)

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


function find_x(model::Yazdani1, cones::AbstractVector{Cone}, pars::Pars)

    ~, i_min = findmin(pars.s)
    x        = cones[i_min].center

    return deepcopy(x)
end


################
# Yazdani 2
################


struct Yazdani2 <: Models end


function find_x(model::Yazdani2, cones::AbstractVector{Cone}, pars::Pars)

    ~, i_min = findmin(pars.h_s)
    x        = cones[i_min].center

    return deepcopy(x)
end


################
# Yazdani 3
################


struct Yazdani3 <: Models end


function find_x(model::Yazdani3, cones::AbstractVector{Cone}, pars::Pars)

    aux      = pars.h_s/maximum(pars.h_s) .+ pars.s/maximum(pars.s)
    ~, i_min = findmin(aux)
    x        = cones[i_min].center

    return deepcopy(x)
end


################
# Utilities
################


maxk_i(a, k) = partialsortperm(a, 1:k, rev=true)
