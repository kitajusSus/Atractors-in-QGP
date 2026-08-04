using Lux
using Zygote
using ForwardDiff
using LinearAlgebra
using Statistics
using StaticArrays
using Printf
using Random
using Optimisers
using AttractorsQGP

"""
    utworz_model_pinn(wymiar_ukryty::Int = 64)
Tworzy bezstanowy model Lux PINN.
"""
function utworz_model_pinn(wymiar_ukryty::Int = 64)
    return Lux.Chain(
        Lux.Dense(3 => wymiar_ukryty, tanh),
        Lux.Dense(wymiar_ukryty => wymiar_ukryty, tanh),
        Lux.Dense(wymiar_ukryty => 2)
    )
end

function przygotuj_dane_treningowe(dane::AbstractMatrix{<:Real})
    n_total_rows = size(dane, 1)
    wystepujace_tau = sort(unique(dane[:, 1]))
    steps_per_evol = length(wystepujace_tau)
    n_evolutions = n_total_rows ÷ steps_per_evol

    X_train = zeros(Float64, 3, n_total_rows)
    Y_train = zeros(Float64, 2, n_total_rows)

    for ev in 1:n_evolutions
        idx_startu = (ev - 1) * steps_per_evol + 1
        T₀ = dane[idx_startu, 2]
        ℛ₀ = dane[idx_startu, 3]

        for step in 1:steps_per_evol
            global_idx = (ev - 1) * steps_per_evol + step
            X_train[1, global_idx] = dane[global_idx, 1]
            X_train[2, global_idx] = T₀
            X_train[3, global_idx] = ℛ₀

            Y_train[1, global_idx] = dane[global_idx, 2]
            Y_train[2, global_idx] = dane[global_idx, 3]
        end
    end
    return X_train, Y_train
end

function calc_physics_residual(model, ps, st, τ::Real, T₀::Real, ℛ₀::Real, hydro_model; eps = 1.0e-4)
    u_plus, _ = model([τ + eps, T₀, ℛ₀], ps, st)
    u_minus, _ = model([τ - eps, T₀, ℛ₀], ps, st)
    dudτ_net = (u_plus .- u_minus) ./ (2 * eps)

    u_pred, _ = model([τ, T₀, ℛ₀], ps, st)
    rhs_out = AttractorsQGP.rhs(u_pred, hydro_model, τ)
    dudτ_physics = [rhs_out[1], rhs_out[2]]

    return sum((dudτ_net .- dudτ_physics) .^ 2)
end

function total_loss_verbose(model, ps, st, X_data, Y_data, hydro_model; λ_phys = 0.5, n_colloc = 100)
    N = size(X_data, 2)

    Y_pred, _ = model(X_data, ps, st)
    loss_data = mean((Y_pred .- Y_data) .^ 2)

    loss_phys = 0.0
    rzeczywiste_colloc = min(n_colloc, N)
    for i in rand(1:N, rzeczywiste_colloc)
        τ, T₀, ℛ₀ = X_data[1, i], X_data[2, i], X_data[3, i]
        loss_phys += calc_physics_residual(model, ps, st, τ, T₀, ℛ₀, hydro_model)
    end
    loss_phys /= rzeczywiste_colloc

    loss_total = loss_data + λ_phys * loss_phys
    return loss_total, loss_data, loss_phys
end

function trenuj_pinn_rygorystycznie!(model, ps, st, X_train, Y_train, hydro_model; epochs = 500, lr = 0.001, λ_phys = 0.5)
    opt = Optimisers.Adam(lr)
    opt_state = Optimisers.setup(opt, ps)

    println("🚀 PINN Trening Lux (Epoki: $epochs, λ_phys: $λ_phys)...")
    println("Epoka | Total Loss | Loss Data | Loss Physics")
    println("-------------------------------------------------")

    for epoch in 1:epochs
        loss_total, grads = Zygote.withgradient(ps) do p
            total_loss_verbose(model, p, st, X_train, Y_train, hydro_model; λ_phys = λ_phys, n_colloc = 200)[1]
        end
        opt_state, ps = Optimisers.update(opt_state, ps, grads[1])

        if epoch % 50 == 0 || epoch == 1
            _, l_d, l_p = total_loss_verbose(model, ps, st, X_train, Y_train, hydro_model; λ_phys = λ_phys, n_colloc = 200)
            @printf("%4d  | %10.6f | %9.6f | %12.6f\n", epoch, loss_total, l_d, l_p)
        end
    end
    println("✅ Trening zakończony!")
    return ps, st
end

function badaj_wymiar_globalnie(model, ps, st, τ_test::Real, X_data_start; n_probek = 50)
    d_eff_list = Float64[]
    N_total = size(X_data_start, 2)

    rzeczywiste_probki = min(n_probek, N_total)
    losowe_indeksy = rand(1:N_total, rzeczywiste_probki)

    for i in losowe_indeksy
        T₀_test = X_data_start[2, i]
        ℛ₀_test = X_data_start[3, i]

        J = ForwardDiff.jacobian(
            x0 -> first(model([τ_test, x0[1], x0[2]], ps, st)),
            [T₀_test, ℛ₀_test]
        )

        J64 = Matrix{Float64}(J)
        C_local = J64 * J64'
        λ = eigvals(Symmetric(C_local))

        if !isempty(λ)
            max_λ = maximum(abs.(λ))
            if max_λ > 0
                λ = λ[λ .> (max_λ * 1.0e-6)] # Względny próg szumu
            end

            if !isempty(λ)
                d_eff = (sum(λ)^2) / sum(λ .^ 2)
                push!(d_eff_list, d_eff)
            end
        end
    end

    return mean(d_eff_list), std(d_eff_list)
end
