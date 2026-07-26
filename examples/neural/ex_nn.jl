using AtractorsQGP
using Lux
using Optimisers
using Zygote
using LinearAlgebra
using Statistics
using ComponentArrays
using Random
using GLMakie

model = HJSWModel()
ics = generate_initial_conditions(model, 100; seed = 5)
sols = generate_trajectories(model, ics, (0.2, 1.5); saveat = 0.01, parallel = :serial)

n_trajs = length(sols)
n_steps = length(sols[1].t)
total_samples = n_trajs * n_steps

X_raw = Matrix{Float64}(undef, total_samples, 4)
Y_raw = Matrix{Float64}(undef, total_samples, 3)

row_idx = 1
for i in 1:n_trajs
    x0 = ics[i]
    for j in 1:n_steps
        X_raw[row_idx, 1] = sols[i].t[j]
        X_raw[row_idx, 2:4] .= x0
        Y_raw[row_idx, :] .= sols[i].u[j]
        row_idx += 1
    end
end

X_norm, mn_in, mx_in = normalize_generic(X_raw)
Y_norm, mn_out, mx_out = normalize_generic(Y_raw)

X_train = Matrix(X_norm')
Y_train = Matrix(Y_norm')

nn = build_generic_network(4, 3; hidden_layers = 3, hidden_size = 64)
rng = Xoshiro(5)
ps_raw, st = Lux.setup(rng, nn)
ps = ComponentArray{Float64}(ps_raw)

opt = Optimisers.Adam(1.0e-3)
opt_state = Optimisers.setup(opt, ps)

println("Training neural network on HJSW data...")
for epoch in 1:2000
    loss_val, (∇ps,) = Zygote.withgradient(ps) do p
        Y_pred, _ = nn(X_train, p, st)
        return mean(abs2, Y_pred .- Y_train)
    end

    global opt_state, ps = Optimisers.update(opt_state, ps, ∇ps)

    if epoch % 100 == 0
        println("Epoch: $epoch | Loss (MSE): $(round(loss_val, sigdigits = 4))")
    end
end

taus_pinn = sols[1].t
x0_ref = ics[1]

dimensions_pinn = Float64[]
for τ in taus_pinn
    J = compute_data_jacobian(nn, ps, st, τ, x0_ref, mn_in, mx_in, mn_out, mx_out)
    σ = svdvals(J)
    d_eff = d_eff_from_singular_values(σ)
    push!(dimensions_pinn, d_eff)
end

wyniki_matrix_neural = hcat(taus_pinn, dimensions_pinn)
fig_neural = plot_pinn_deff_evolution(wyniki_matrix_neural)

out_file = joinpath(@__DIR__, "hjsw_neural_deff_evolution.png")
save(out_file, fig_neural)
println("Saved plot: ", out_file)
