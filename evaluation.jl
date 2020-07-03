function compute_performance!(T_surv::Int, δ::Real, model::Models, cones::AbstractVector{Cone}, pars::Pars)

    surv = zeros(T_surv)

    x = find_x(model, cones, pars; δ=δ)
    for t in 1:T_surv

        update_cones!(cones, pars)

        if t==1 || cones_value(x, cones) < δ
            x = find_x(model, cones, pars; δ=δ)

            surv[t] = 0
        else
            surv[t] = surv[t-1] + 1
        end
    end

    return mean(surv)
end


function compute_performance(T_surv::Int, δ::Real, model::AbstractVector{Models}, cones::AbstractVector{Cone}, pars::Pars)

    return compute_performance!.(T_surv, δ, models, [deepcopy(cones)], [pars])
end


function compute_performance_enhanced!(T_surv::Int, δ::Real, model::Models, cones::AbstractVector{Cone}, pars::Pars)

    surv = zeros(T_surv)

    x = find_x(model, cones, pars; δ=δ)
    for t in 1:T_surv

        update_cones!(cones, pars)

        if t==1 || cones_value(x, cones) < δ
            x = find_x(model, cones, pars; δ=δ)

            surv[t] = 0
        else
            surv[t] = surv[t-1] + 1
        end
    end

    ii_surv  = findall(surv .== 0)
    surv_mod = surv[[ii_surv[2:end].-1; length(surv)]]

    return [mean(surv), mean(surv_mod)]
end


function compute_performance_enhanced(T_surv::Int, δ::Real, model::AbstractVector{Models}, cones::AbstractVector{Cone}, pars::Pars)

    return hcat(compute_performance_enhanced!.(T_surv, δ, models, [deepcopy(cones)], [pars])...)
end


function compute_metrics!(x::AbstractMatrix, T_surv::Int, T_aver::Int, δ::Real, cones::AbstractVector{Cone}, pars::Pars)

    n      = size(x)[2]
    f_surv = zeros(n)
    f_aver = zeros(n)
    ii     = trues(n)

    update_cones!(cones, pars) # Update cones to start from the next time instant

    for t in 1:max(T_aver, T_surv)

        f = cones_value(x, cones)
        if t <= T_surv
            ii     .*= (f .>= δ)
            f_surv .+= ii
        end
        if t <= T_aver
            f_aver .+= f
        end

        update_cones!(cones, pars)

        t > T_aver && sum(ii) == 0 && break
    end
    return f_surv, f_aver./T_aver
end


function compute_metrics(n_try::Int, x::AbstractMatrix, T_surv::Int, T_aver::Int, δ::Real, cones::AbstractVector{Cone}, pars::Pars)

    f_all  = [compute_metrics!(x, T_surv, T_aver, δ, deepcopy(cones), pars) for i in 1:n_try]
    f_surv = hcat([f_try[1] for f_try in f_all]...)
    f_aver = hcat([f_try[2] for f_try in qwe]...)
    return mean(f_surv, dims=2), mean(f_aver, dims=2)
end


function compute_region_attraction(n::Int, cones::AbstractVector{Cone}, pars::Pars)

    x      = random_points(n, pars)
    i_cone = cones_max_i(x, cones)
    counts = [sum(i_cone .== i) for i in 1:pars.m]
end
