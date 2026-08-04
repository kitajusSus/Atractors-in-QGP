using AttractorsQGP
using GLMakie
using Statistics
using Printf

include(joinpath(@__DIR__, "../local_pca/lpca.jl"))
include(joinpath(@__DIR__, "../local_pca/ex_lpca.jl"))
include(joinpath(@__DIR__, "functions_pinn.jl"))
include(joinpath(@__DIR__, "../lle_examples/swissroll.jl"))
include(joinpath(@__DIR__, "../lle_examples/helix.jl"))
include(joinpath(@__DIR__, "../lle_examples/scurve.jl"))

function przeprowadz_testy_geometryczne()
    println("\n--- ETAP 1: L-PCA (Swiss Roll) ---")
    X_swiss, _ = swiss_roll(1500, noise = 0.0)

    for k in [10, 30, 80, 150]
        sredni_wymiar = mean(dims(X_swiss, k = k, tol = 0.05))
        @printf("K = %3d | Średni wymiar: %.2f\n", k, sredni_wymiar)
    end
    return
end

function badaj_qgp_klasycznie(dane)
    println("\n--- ETAP 2: L-PCA (QGP) ---")
    fig = main_local_pca(dane; n_slices = 20, tablica_k = [10, 20, 40], feature_cols = [2, 3])

    save(joinpath(@__DIR__, "lpca_qgp_comparison.png"), fig)
    return println("Zapisano: lpca_qgp_comparison.png")
end

function badaj_qgp_pinn(dane)
    println("\n--- ETAP 3: Rygorystyczny PINN (QGP) ---")
    X_train, Y_train = przygotuj_dane_treningowe(dane)

    fizyka = MISModel()
    model_pinn = utworz_model_pinn(64)

    trenuj_pinn_rygorystycznie!(model_pinn, X_train, Y_train, fizyka; epochs = 100, lr = 0.001, λ_phys = 1)

    return analizuj_zanik_wymiaru(model_pinn, X_train)
end

function analizuj_zanik_wymiaru(model, X_train)

    println("\n" * "="^60)
    println("🔍 GĘSTY SKAN: SZUKANIE CZASU KOLAPSU WYMIARU (ATRAKTORA)")
    println("="^60)

    # 1. Definiujemy gęstą siatkę czasu tau
    taus = collect(range(0.2, 1.5, length = 60))
    srednie_d_eff = Float64[]
    odchylenia_d_eff = Float64[]

    print("Skanowanie czasu τ w toku (60 punktów)... ")
    for tau in taus
        # Korzystamy z funkcji z functions_pinn.jl
        sr, std_err = badaj_wymiar_globalnie(model, tau, X_train, n_probek = 80)
        push!(srednie_d_eff, sr)
        push!(odchylenia_d_eff, std_err)
        print(".") # Pasek postępu
    end
    println(" Zakończono!")

    # 2. Matematyczne szukanie czasu kolapsu (tau_c)
    prog_atraktora = 1.05 # Definiujemy, że d_eff < 1.05 to wymiar 1D
    idx_kolapsu = findfirst(x -> x < prog_atraktora, srednie_d_eff)

    # Note: tau_c is set here for visualization compatibility
    tau_c = isnothing(idx_kolapsu) ? NaN : taus[idx_kolapsu]

    fig = Figure(size = (800, 500))
    ax = Axis(
        fig[1, 1],
        title = "Ewolucja Efektywnego Wymiaru Atraktora (PINN)",
        xlabel = "Czas własny τ [fm/c]",
        ylabel = "Efektywny Wymiar d_eff (Participation Ratio)"
    )

    band!(
        ax, taus,
        srednie_d_eff .- odchylenia_d_eff,
        srednie_d_eff .+ odchylenia_d_eff,
        color = (:blue, 0.2)
    )

    lines!(ax, taus, srednie_d_eff, color = :blue, linewidth = 3, label = "Wymiar d_eff(τ)")

    # Linie pomocnicze (Idealne 2D, Idealne 1D i Czas Kolapsu)
    hlines!(ax, [2.0], color = :gray, linestyle = :dash, label = "Wymiar równowagi 2D")
    hlines!(ax, [1.0], color = :red, linestyle = :dash, label = "Atraktor 1D")

    if !isnan(tau_c)
        vlines!(ax, [tau_c], color = :green, linestyle = :dashdot, label = "Czas kolapsu τ_c ≈ $(round(tau_c, digits = 2))")
    end

    axislegend(ax)

    # Zapis
    save(joinpath(@__DIR__, "zanik_wymiaru_pinn.png"), fig)
    println("  'zanik_wymiaru_pinn.png'")

    return tau_c
end
function uruchom_pelna_analize()
    sciezka = "datasets/ncbj_xin_2020_5000.h5"

    dane = load_dataset(sciezka)

    przeprowadz_testy_geometryczne()
    badaj_qgp_klasycznie(dane)
    return badaj_qgp_pinn(dane)
end

# Uruchomienie całości
uruchom_pelna_analize()
