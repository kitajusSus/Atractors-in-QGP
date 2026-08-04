using LinearAlgebra
using Statistics
using GLMakie
using AttractorsQGP

GLMakie.activate!()

include("lle_examples/helix.jl")
include("lle_examples/swissroll.jl")
include("lle_examples/scurve.jl")
include("lle_examples/2d_examples_data.jl")
include("local_pca/lpca.jl")

function przygotuj_zbiory_danych(N::Int)
    zestawy = [
        ("Płaski Dywan", flat_carpet_dane),
        ("Swiss Roll", swissroll_dane),
        ("Krzywa S", scurve_dane),
        ("Sinus", sinus_dane),
        ("Hiperbola", hiperbola_dane),
    ]

    spakowane_dane = []

    for (nazwa, funkcja) in zestawy
        X, labels = funkcja(N)
        push!(spakowane_dane, (nazwa = nazwa, X = X, labels = labels))
    end

    return spakowane_dane
end

function symulacja_dims(N, K)
    spakowane_dane = przygotuj_zbiory_danych(N)
    ℵ = Vector{Vector{Float64}}()
    nazwy = String[]

    for dane in spakowane_dane
        W = ncbj4_lle_basic(dane.X, K)
        M = (I - W)' * (I - W)
        F = eigen(Symmetric(M))

        println(dane.nazwa)
        println("----------")
        push!(ℵ, F.values)
        push!(nazwy, dane.nazwa)
    end

    for dane in spakowane_dane
        dd = dims(dane.X, k = K)
        meandim = mean(dd)
        println("$(dane.nazwa) - dla K = $(K) Średnia wymiarowość: $meandim")
    end

    return ℵ, nazwy
end

function compare_eigen(ℵ::Vector{Vector{Float64}}, num_evals::Int = 15)
    ℵ_skalowane = Vector{Vector{Float64}}()

    for λ in ℵ
        λ_sorted = sort(real.(λ))
        limit = min(num_evals, length(λ_sorted))
        λ_sub = copy(λ_sorted[1:limit])

        if length(λ_sub) >= 2 && abs(λ_sub[2]) > 1.0e-15
            λ_sub[2:end] .= λ_sub[2:end] ./ λ_sub[2]
        end

        push!(ℵ_skalowane, λ_sub)
    end

    return ℵ_skalowane
end

function plot_compare_dims(ℵ::Vector{Vector{Float64}}, nazwy::Vector{String})
    set_publication_theme_large()

    num_evals = 5
    ℵ_skalowane = compare_eigen(ℵ, num_evals)

    fig = Figure(size = (1200, 800))
    ax = Axis(
        fig[1, 1],
        title = "Porównanie znormalizowanego widma LLE",
        xlabel = "Indeks wartości własnej",
        ylabel = "Wartość znormalizowana (skala log)",
        xticks = 1:num_evals
    )

    for i in 1:length(ℵ_skalowane)
        x_vals = 1:length(ℵ_skalowane[i])
        y_vals = max.(ℵ_skalowane[i], 1.0e-15)

        lines!(ax, x_vals, y_vals, label = nazwy[i])
        scatter!(ax, x_vals, y_vals)
    end

    axislegend(ax, position = :rb)

    return fig
end

function analiza_skokow_K(X, K_range, max_d = 6)
    ratios = zeros(Float64, length(K_range), max_d)

    for (i, K) in enumerate(K_range)
        W = ncbj4_lle_basic(X, K)
        M = (I - W)' * (I - W)
        F = eigen(Symmetric(M))

        λ = sort(real.(F.values))

        for d in 1:max_d
            val_down = max(λ[d + 1], 1.0e-15)
            val_up = max(λ[d + 2], 1.0e-15)
            ratios[i, d] = val_up / val_down
        end
    end

    mean_ratios = vec(mean(ratios, dims = 1))
    var_ratios = vec(var(ratios, dims = 1))

    return ratios, mean_ratios, var_ratios
end

function plot_analiza_skokow(nazwa, K_range, ratios, mean_ratios, var_ratios, max_d)
    set_publication_theme_large()

    fig = Figure(size = (1600, 800))

    ax1 = Axis(
        fig[1, 1],
        title = "[$nazwa] Stosunek wartości (λ_d+2 / λ_d+1) od K",
        xlabel = "Liczba sąsiadów K",
        ylabel = "Stosunek (skala log)",
        yscale = log10
    )

    for d in 1:max_d
        lines!(ax1, K_range, ratios[:, d], label = "d = $d")
        scatter!(ax1, K_range, ratios[:, d], markersize = 10)
    end
    axislegend(ax1, position = :rt)

    ax2 = Axis(
        fig[1, 2],
        title = "[$nazwa] Statystyka skoków dla K ∈ [$(first(K_range)), $(last(K_range))]",
        xlabel = "Kandydat na wymiar (d)",
        ylabel = "Wartość (skala log)",
        yscale = log10,
        xticks = 1:max_d
    )

    barplot!(
        ax2, 1:max_d, max.(mean_ratios, 1.0e-10),
        color = (:dodgerblue, 0.7), label = "Średni skok (Mean)"
    )

    scatterlines!(
        ax2, 1:max_d, max.(var_ratios, 1.0e-10),
        color = :crimson, marker = :diamond, markersize = 20, linewidth = 4, label = "Wariancja (Var)"
    )

    axislegend(ax2, position = :rt)

    return fig
end

function przeprowadz_pelna_analize(N = 1500, K_range = 4:25, max_d = 6)
    spakowane_dane = przygotuj_zbiory_danych(N)
    wyniki_koncowe = Dict{String, Any}()

    for dane in spakowane_dane
        println("Przetwarzanie zbioru: $(dane.nazwa)")

        ratios, mean_ratios, var_ratios = analiza_skokow_K(dane.X, K_range, max_d)
        fig = plot_analiza_skokow(dane.nazwa, K_range, ratios, mean_ratios, var_ratios, max_d)

        display(fig)

        wyniki_koncowe[dane.nazwa] = (
            X = dane.X,
            labels = dane.labels,
            ratios = ratios,
            mean_ratios = mean_ratios,
            var_ratios = var_ratios,
            fig = fig,
        )
    end

    return wyniki_koncowe
end


function run_analiza()

    paczka = przeprowadz_pelna_analize(2000, 5:30, 6)
    for i in 1:30
        for (nazwa, dane) in paczka
            println("Zbiór: $(nazwa)")
            println("Średnie skoki: $(dane.mean_ratios)")
            println("Wariancje skoków: $(dane.var_ratios)")
            println("--------------------------------------------------")
        end
        aleph, nazwy = symulacja_dims(2000, i)
        fig_widmo = plot_compare_dims(aleph, nazwy)
        display(fig_widmo)
    end

    return paczka, fig_widmo
end


# if abspath(PROGRAM_FILE) == @__FILE__
#     paczka = przeprowadz_pelna_analize(2000, 5:30, 6)
#
#     aleph, nazwy = symulacja_dims(2000, 12)
#     fig_widmo = plot_compare_dims(aleph, nazwy)
#     display(fig_widmo)
#
#     readline()
# end
