include("pars.jl")
include("cone.jl")
include("models.jl")
include("evaluation.jl")
include("utilities.jl")

n_eval = 2500

T_surv = 100
T_aver = 4
δ      = 40

d = 2
#m = 10
m = 14

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

models     = [TMO(); Random(); RandomAbove(); Robust1(); RoA1(n_eval); Yazdani1(); Yazdani2(); Yazdani3(); Yazdani4(); Yazdani5()]
n_models = length(models)
#models     = [TMO(); Robust1(); Yazdani1(); Yazdani2(); Yazdani3(); Yazdani4(); Yazdani5()]


"""
#n_try = 1000
n_try = 1000

surv1 = zeros(length(models), n_try)
surv2 = zeros(length(models), n_try)
for i in 1:n_try

    mod(i, 100) == 0 && println("Iteration " * string(i))

    pars_fixed = Pars(d, m, x_min, x_max, h_min, h_max, w_min, w_max, h_init, w_init, s_min, s_max, h_s_min, h_s_max, w_s_min, w_s_max, generator)
    cones      = init_cones(pars_fixed)

    surv_all   = compute_performance_enhanced(T, δ, models, cones, pars_fixed)
    surv1[:,i] = surv_all[1,:]
    surv2[:,i] = surv_all[2,:]
end

[mean(surv1, dims=2) mean(surv2, dims=2)]
"""


"""
δ_n = 25
surv_red1 = hcat([mean(surv1[:,n_beg+1:n_beg+δ_n], dims=2) for n_beg in 0:δ_n:(n_try-δ_n)]...)
δ_n = 100
surv_red2 = hcat([mean(surv1[:,n_beg+1:n_beg+δ_n], dims=2) for n_beg in 0:δ_n:(n_try-δ_n)]...)



histogram(surv_red1[1,:], normalize=true, bins=5)
histogram!(surv_red2[1,:], normalize=true, bins=5)
"""








n_try = 1000


data = zeros(n_models+2, m, n_try)
for i in 1:n_try

    pars    = Pars(d, m, x_min, x_max, h_min, h_max, w_min, w_max, "uniform", "uniform", s_min, s_max, h_s_min, h_s_max, w_s_min, w_s_max, generator)
    cones   = init_cones(pars)
    centers = center(cones)

    data[1,:,i], data[2,:,i] = compute_metrics!(centers, T_surv, T_aver, δ, deepcopy(cones), pars)
    data[3:end,:,i]          = hcat(scores.(models, [cones], [pars]; δ=δ)...)'
end

data_mean = dropdims(mean(data, dims=3), dims=3)

# preved na ranky


using StatsBase

n_data = n_models+2
corrs  = zeros(n_data, n_data, n_try)
for i in 1:n_data
    for j in 1:n_data
        for k in 1:n_try
            corrs[i,j,k] = cor(tiedrank(data[i,:,k], rev=true), tiedrank(data[j,:,k], rev=true))
        end
    end
end




corrs_mean = zeros(n_data, n_data)
for i in 1:n_data
    for j in 1:n_data
        corrs_mean[i,j] = mean(corrs[i,j,.!isnan.(corrs[i,j,:])])
    end
end




round.(corrs_mean[:,1:2]; digits=2)




foijioef


1.0    0.77
0.77   1.0
0.51   0.42
-0.01  -0.0
0.01   0.01
0.58   0.51
0.28   0.31
0.5    0.41
0.25   0.31
0.05  -0.0
0.2    0.2
0.42   0.37



1.0   0.77
0.77  1.0
0.51  0.43
0.01  0.01
0.02  0.01
0.58  0.5
0.29  0.31
0.51  0.41
0.25  0.32
0.06  0.01
0.2   0.21
0.43  0.38





a = 1
