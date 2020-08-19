using Printf
using BSON: @save

include("pars.jl")
include("cone.jl")
include("models.jl")
include("evaluation.jl")
include("utilities.jl")
include("plots.jl")

folder_name = "Results_Final"

δ_all  = collect(20:2:50)
n_try  = 100000
n_δ    = length(δ_all)

T_surv = 100
T_aver = 4

δ = 40
d = 5
m = 10

x_min     = -50.
x_max     = 50.
h_min     = 30.
h_max     = 70.
w_min     = 1.
w_max     = 12.
s_min     = 0.5
s_max     = 3.
h_s_min   = 1.
h_s_max   = 15.
w_s_min   = 0.1
w_s_max   = 1.5
generator = "Normal"

models   = [RoA1(2500); TMO(); Robust1(); Robust2(); Yazdani1(); Yazdani2(); Yazdani3(); Yazdani4(); RandomAbove(); Random()]
n_models = length(models)
n_data   = n_models+2

pars = Pars(d, m, x_min, x_max, h_min, h_max, w_min, w_max, "uniform", "uniform", s_min, s_max, h_s_min, h_s_max, w_s_min, w_s_max, generator)
create_directory(folder_name)
for (i_δ, δ) in enumerate(δ_all)

    @printf("Running experiment %3d out of %d\n", i_δ, n_δ)

    data = zeros(n_data, m, n_try)
    for i in 1:n_try

        pars    = Pars(d, m, x_min, x_max, h_min, h_max, w_min, w_max, "uniform", "uniform", s_min, s_max, h_s_min, h_s_max, w_s_min, w_s_max, generator)
        cones   = init_cones(pars)
        centers = center(cones)

        data[1,:,i], data[2,:,i] = compute_metrics!(centers, T_surv, T_aver, δ, deepcopy(cones), pars)
        data[3:end,:,i]          = hcat(scores.(models, [cones], [pars]; δ=δ)...)'
    end

    metrics = zeros(2, n_models, n_try)
    for i in 1:n_try, k in 3:n_data

        i_max            = findmax(data[k,:,i])[2]
        metrics[:,k-2,i] = data[1:2,i_max,i]
    end

    pars      = Pars(d, m, x_min, x_max, h_min, h_max, w_min, w_max, "uniform", "uniform", s_min, s_max, h_s_min, h_s_max, w_s_min, w_s_max, generator)
    file_name = get_file_name(folder_name, δ, m)
    @save file_name data metrics pars models T_surv T_aver
end
