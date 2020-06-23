random_points(n::Int, pars::Pars) = pars.x_min .+ (pars.x_max .- pars.x_min).*rand(pars.d, n)
