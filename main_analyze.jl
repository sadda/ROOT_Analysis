using Plots
using StatsBase
using DelimitedFiles
using KernelDensity
using BSON: @load

include("pars.jl")
include("cone.jl")
include("models.jl")
include("evaluation.jl")
include("utilities.jl")
include("plots.jl")

macro R_str(s)
    s
end

δ = 40
m = 10
pars_all = collect(20:2:50)

header_metrics  = [R"$\Fsurv$", R"$\Faver$"]
header_methods  = ["TMO" "Random" "RandomAbove" "Robust" "RoA" "Yazdani1" "Yazdani2" "Yazdani3" "Yazdani4" "Yazdani5"]
header_notation = [R"$d$", R"$m$", R"$\delta$", R"$x_{\rm min}$", R"$x_{\rm max}$", R"$h_{\rm min}$", R"$h_{\rm max}$", R"$w_{\rm min}$", R"$w_{\rm max}$", R"$s_{\rm min}$", R"$s_{\rm max}$", R"$\sigma_{h,\rm min}$", R"$\sigma_{h,\rm max}$", R"$\sigma_{w,\rm min}$", R"$\sigma_{w,\rm max}$"]
ii_red          = [1, 4, 5, 6, 7, 10]


## Prepare data


@load get_file_name(δ, m) data metrics pars models T_surv T_aver

n_pars  = length(pars_all)
n_model = size(metrics)[2]
n_try   = size(metrics)[3]
f_surv  = zeros(n_model, n_pars)
f_aver  = zeros(n_model, n_pars)
corrs   = zeros(n_model, 2, n_pars, n_try)
for (i, δ) in enumerate(pars_all)

    @printf("Computing correlations %3d out of %3d\n", i, n_pars)
    @load get_file_name(δ, m) data metrics pars models T_surv T_aver

    f_surv[:,i] = dropdims(mean(metrics, dims=3), dims=3)[1,:]
    f_aver[:,i] = dropdims(mean(metrics, dims=3), dims=3)[2,:]

    for j in 1:n_model, k in 1:2
        for l in 1:n_try
            corrs[j,k,i,l] = cor(tiedrank(data[j+2,:,l], rev=true), tiedrank(data[k,:,l], rev=true))
        end
    end
end

corrs_mean = zeros(n_model, 2)
for j in 1:n_model, k in 1:2
    corrs_aux       = reshape(corrs[j,k,:,:], n_pars*n_try, 1)
    corrs_mean[j,k] = mean(corrs_aux[.!isnan.(corrs_aux)])
end


## Compute correlations


bcorrs      = zeros(n_model, 2, n_try)
bcorrs_mean = zeros(n_model, 2)
@time for i in 1:n_model, j in 1:2
    for k in 1:n_try
        bcorrs[i,j,k] = cor(tiedrank(data[i+2,:,k], rev=true), tiedrank(data[j,:,k], rev=true))
    end
    bcorrs_mean[i,j] = mean(bcorrs[i,j,.!isnan.(bcorrs[i,j,:])])
end




acorrs      = zeros(n_model, 2, n_try)
acorrs_mean = zeros(n_model, 2)
@time for i in 1:n_model, j in 1:2
    acorrs_mean[i,j] = mean(corrs[i,j,11,.!isnan.(corrs[i,j,11,:])])
end

## Make plots


p1 = plot(pars_all, f_surv[ii_red,:]', label=header_methods[ii_red'])

surv_rank_red = hcat([tiedrank(f_surv[ii_red,k], rev=true) for k in 1:n_pars]...)
p2 = plot(pars_all, surv_rank_red', legend=:none)

p3 = plot(pars_all, f_aver[ii_red,:]', legend=:none)

aver_rank_red = hcat([tiedrank(f_aver[ii_red,k], rev=true) for k in 1:n_pars]...)
p4 = plot(pars_all, aver_rank_red', legend=:none)

p5 = plot(p1, p2, p3, p4)
savefig(p5, "F_all.png")


## Plot histograms for different repetition numbers


compute_rolling_mean(x, δ_n) = hcat([mean(x[:,n_beg+1:n_beg+δ_n], dims=2) for n_beg in 0:δ_n:(size(x)[2]-δ_n)]...)


function compute_density(x::AbstractArray)

    n_x    = length(x)
    dens   = [kde(x[i]) for i in 1:n_x]
    f_min  = minimum([minimum(dens[i].x) for i in 1:n_x])
    f_max  = maximum([maximum(dens[i].x) for i in 1:n_x])
    x_step = (maximum(dens[1].x) - minimum(dens[1].x)) / length(dens[1].x)

    res = Array{Array{Float64}}(undef, n_x)
    for i in 1:n_x
        res[i] = hcat([f_min-x_step; collect(dens[i].x); f_max-x_step], [0; dens[i].density; 0])
    end

    return res
end


n_method = 10
δ_all    = [400 1000 4000]
surv_red = [compute_rolling_mean(metrics[1,:,:], δ_n)[n_method,:] for δ_n in δ_all]
dens     = compute_density(surv_red)

p6 = plot_repetition(dens)
savefig(p6, "Hist.png")


## Plot region at attraction


d = 2
m = 3
pars_mod = Pars(d, m, pars.x_min[1], pars.x_max[1], pars.h_min, pars.h_max, pars.w_min, pars.w_max, "uniform", "uniform", 1., 1., 1., 1., 1., 1., "Normal")
cones    = init_cones(pars_mod)
cones[1].center = [-30, -30]
cones[1].height = 50
cones[1].width  = pars.w_min
cones[2].center = [30, -30]
cones[2].height = 70
cones[2].width  = pars.w_max
cones[3].center = [0, 30]
cones[3].height = 60
cones[3].width  = 0.5*(pars.w_min+pars.w_max)

p7 = plot_contours(cones, pars_mod; min_value=-Inf, plot_boundary=true)
savefig(p7, "Contours.png")


## Save results to csv


writedlm("surv1.csv",  hcat(pars_all, f_surv[ii_red,:]'), ',')
writedlm("surv2.csv",  hcat(pars_all, surv_rank_red'), ',')
writedlm("aver1.csv",  hcat(pars_all, f_aver[ii_red,:]'), ',')
writedlm("aver2.csv",  hcat(pars_all, aver_rank_red'), ',')
writedlm("hist.csv",  hcat([dens[i] for i in 1:length(dens)]...), ',')


## Print Latex tables


caption  = "Correlations between models and metrics"
label    = "table:corrs"

table_to_tex(corrs_mean, "%1.1p"; header_l=header_methods, header_t=header_metrics, caption=caption, label=label, alignment=alignment_lr(header_metrics))


"""
caption  = "Survival time and averaged objective metrics"
label    = "table:metrics"

metrics_red = dropdims(mean(metrics, dims=3), dims=3)'
table_to_tex(metrics_red, "%1.1f"; header_l=header_methods, header_t=header_metrics, caption=caption, label=label, alignment=alignment_lr(header_metrics))
"""


caption  = "Default hyperparameters"
label    = "table:pars"


s_min     = 0.5
s_max     = 3.
h_s_min   = 1.
h_s_max   = 15.
w_s_min   = 0.1
w_s_max   = 1.5
data_notation = [d, m, δ, pars.x_min[1], pars.x_max[1], pars.h_min, pars.h_max, pars.w_min, pars.w_max, s_min, s_max, h_s_min, h_s_max, w_s_min, w_s_max]
table_to_tex(data_notation, "%1.1f"; header_l=header_notation, caption=caption, label=label, alignment=alignment_lr(ones(1)))
