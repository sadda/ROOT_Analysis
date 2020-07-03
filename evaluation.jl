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
    f_aver = zeros(n)
    f_surv = zeros(n)
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
