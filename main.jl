include("pars.jl")
include("cone.jl")
include("models.jl")
include("evaluation.jl")
include("utilities.jl")
include("plots.jl")

n_eval = 2500

T_surv = 100
T_aver = 4
δ      = 40

d = 2
m = 20

x_min     = -50.
x_max     = 50.
h_min     = 30.
h_max     = 70.
w_min     = 1.
w_max     = 12.
s_min     = 0.5
s_max     = 3.
h_init    = 50.
w_init    = 6.
h_s_min   = 1.
h_s_max   = 15.
w_s_min   = 0.1
w_s_max   = 1.5
generator = "Normal"

models   = [TMO(); Random(); RandomAbove(); Robust1(); RoA1(n_eval); Yazdani1(); Yazdani2(); Yazdani3(); Yazdani4(); Yazdani5()]
n_models = length(models)
n_data   = n_models+2


## Compute survival time and averaged objective


n_try = 10000

data = zeros(n_data, m, n_try)
for i in 1:n_try

    pars    = Pars(d, m, x_min, x_max, h_min, h_max, w_min, w_max, "uniform", "uniform", s_min, s_max, h_s_min, h_s_max, w_s_min, w_s_max, generator)
    cones   = init_cones(pars)
    centers = center(cones)

    data[1,:,i], data[2,:,i] = compute_metrics!(centers, T_surv, T_aver, δ, deepcopy(cones), pars)
    data[3:end,:,i]          = hcat(scores.(models, [cones], [pars]; δ=δ)...)'
end


## Compute performance of models


metrics = zeros(2, n_models, n_try)
for i in 1:n_try, k in 3:n_data

    i_max            = findmax(data[k,:,i])[2]
    metrics[:,k-2,i] = data[1:2,i_max,i]
end


## Compute ranks and correlations between methods


using StatsBase

corrs      = zeros(n_data, n_data, n_try)
corrs_mean = zeros(n_data, n_data)
for i in 1:n_data, j in 1:n_data
    for k in 1:n_try
        corrs[i,j,k] = cor(tiedrank(data[i,:,k], rev=true), tiedrank(data[j,:,k], rev=true))
    end
    corrs_mean[i,j] = mean(corrs[i,j,.!isnan.(corrs[i,j,:])])
end


## Plot histograms for different repetition numbers


compute_rolling_mean(x, δ_n) = hcat([mean(x[:,n_beg+1:n_beg+δ_n], dims=2) for n_beg in 0:δ_n:(size(x)[2]-δ_n)]...)

surv_red1 = compute_rolling_mean(metrics[1,:,:], 250)
surv_red2 = compute_rolling_mean(metrics[1,:,:], 1000)

n_method = 10
δ_all    = [250 1000]
surv_red = [compute_rolling_mean(metrics[1,:,:], δ_n)[n_method,:] for δ_n in δ_all]

plot_repetition(surv_red)


## Print correlations in a Latex table


macro R_str(s)
    s
end

caption  = "Correlations between models (rows) and metrics (columns)."
label    = "fig:corrs"
header_t = [R"$\Fsurv$", R"$\Faver$"]
header_l = ["TMO", "Random", "RandomAbove", "Robust", "RoA", "Yazdani1", "Yazdani2", "Yazdani3", "Yazdani4", "Yazdani5"]

table_to_tex(corrs_mean[3:end,1:2], "%1.1p"; header_l=header_l, header_t=header_t, caption=caption, label=label)
















a = 1
