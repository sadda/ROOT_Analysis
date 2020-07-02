include("pars.jl")
include("cone.jl")
include("models.jl")
include("evaluation.jl")
include("utilities.jl")

n_eval = 2500

T = 100
δ = 40

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

pars_fixed = Pars(d, m, x_min, x_max, h_min, h_max, w_min, w_max, h_init, w_init, s_min, s_max, h_s_min, h_s_max, w_s_min, w_s_max, generator)
pars_uni   = Pars(d, m, x_min, x_max, h_min, h_max, w_min, w_max, "uniform", "uniform", 1., 1., h_s_min, h_s_max, w_s_min, w_s_max, generator)
models     = [TMO(); Robust1(n_eval); RoA(); Yazdani1(); Yazdani2(); Yazdani3(); Yazdani4(); Yazdani5()]

n_try = 1000
#n_try = 10
res   = zeros(n_try)

surv = zeros(length(models), n_try)

for i in 1:n_try

    mod(i, 100) == 0 && println("Iteration " * string(i))

    cones     = init_cones(pars_fixed)
    surv[:,i] = compute_performance(T, δ, models, cones, pars_fixed)
end

println(mean(surv, dims=2))






"""
cones   = init_cones(pars_uni)
centers = center(cones)
heights = height(cones)
widths  = width(cones)

f_surv, f_aver  = compute_metrics(50, centers, T, 10, δ, cones, pars_uni)
counts          = compute_region_attraction(10000, cones, pars_uni)
aux             = heights/maximum(heights) .- widths/maximum(widths)

data = hcat(counts/sum(counts), f_surv, f_aver, aux, heights, widths)
data = sortslices(data, dims=1; rev=true)

data
"""
