
include("lib.jl")
using .modHydroSim
using Plots
using Statistics
using LinearAlgebra

module PCAWorkflow

using ..modHydroSim
using Plots
using Statistics
using LinearAlgebra

export SimSettings, run_simulation, analyze_pca_over_time, plot_all_pca_results, plot_pca_at_time

function pca_math(X::Matrix{Float64}, n_components::Int)
    # Pobranie wymiarów macierzy wejściowej
    n_samples, n_features = size(X) # n = n_samples, p = n_features
# 
    if n_components > n_features
        error("Liczba komponentów (k) nie może być większa niż liczba cech (p).")
    end

    # --------------------------------------------------------------------------
    # Krok 1: Centrowanie danych
    #
    # Cel: Przesunięcie danych tak, aby każda cecha miała średnią równą 0.
    # Wzór: Dla każdej kolumny (cechy) X_j, obliczana jest jej średnia miu_j.
    #       Następnie dla każdej obserwacji x_ij odejmuje się tę średnią:
    #       x'_ij = x_ij - miu_j
    # W notacji macierzowej: X_c = X - μ (gdzie μ to wektor średnich powielony dla każdego wiersza)
    # --------------------------------------------------------------------------
    mean_vector = mean(X, dims=1) # Wektor średnich dla każdego wymiaru  (1 x p)
    X_centered = X .- mean_vector

    # -------------------------------------
    # Krok 2: Obliczenie macierzy kowariancji
    #
    # Cel: Stworzenie macierzy (p x p), która opisuje wariancję i współzmienność cech.
    # Wzór: Macierz kowariancji C jest zdefiniowana jako:
    #       C = (1 / (n-1)) * X_c' * X_c
    #       gdzie X_c' to transpozycja macierzy scentrowanych danych.
    #       Element c_ij macierzy C to kowariancja między cechą i a cechą j.
    #       cov(X_i, X_j) = E[(X_i - E[X_i])(X_j - E[X_j])]
    # --------------------------------------------------------------------------
    cov_matrix = cov(X_centered) # Rozmiar (p x p)

    # --------------------------------------------------------------------------
    # Krok 3: Obliczenie wartości własnych i wektorów własnych
    #
    # Cel: Znalezienie głównych osi wariancji w danych.
    # Wzór: Rozwiązujemy równanie własne dla macierzy kowariancji C:
    #       C * v = λ * v
    #       gdzie:
    #       - v to wektor własny macierzy C (kierunek osi, czyli główna składowa).
    #       - λ to wartość własna odpowiadająca wektorowi v (mówi o tym,
    #         jak dużo wariancji jest w kierunku wektora v).
    # To jest równoważne znalezieniu rozwiązań równania charakterystycznego:
    #       det(C - λI) = 0
    eigen_result = eigen(cov_matrix)
    eigenvalues = eigen_result.values    # Wektor λ (lambda)
    eigenvectors = eigen_result.vectors  # Macierz, której kolumnami są wektory v

    # --------------------------------------------------------------------------
    # Krok 4: Sortowanie wektorów własnych i wybór głównych składowych
    #
    # Cel: Uporządkowanie osi wariancji od najważniejszej do najmniej ważnej
    #      i wybranie k pierwszych z nich.
    # Proces: Sortujemy wartości własne λ w porządku malejącym:
    #         λ_1 ≥ λ_2 ≥ ... ≥ λ_p
    #         Następnie, w tej samej kolejności, sortujemy odpowiadające im
    #         wektory własne v.
    #         Wybieramy k pierwszych wektorów własnych (v_1, v_2, ..., v_k),
    #         aby utworzyć macierz transformacji W.
    # --------------------------------------------------------------------------
    sorted_indices = sortperm(eigenvalues, rev=true)
    sorted_eigenvalues = eigenvalues[sorted_indices]
    sorted_eigenvectors = eigenvectors[:, sorted_indices]

    # Macierz transformacji W (projection_matrix) składa się z k pierwszych wektorów własnych
    # W = [v_1 | v_2 | ... | v_k]
    projection_matrix = sorted_eigenvectors[:, 1:n_components] # Rozmiar (p x k)

    # --------------------------------------------------------------------------
    # Krok 5: Transformacja danych do nowej przestrzeni
    #
    # Cel: Rzutowanie oryginalnych danych na nową podprzestrzeń zdefiniowaną
    #      przez wybrane główne składowe.
    # Wzór: Nowa macierz danych Y jest obliczana jako iloczyn macierzy
    #       scentrowanych danych X_c i macierzy transformacji W.
    #       Y = X_c * W
    #       Wynikowa macierz Y ma wymiary (n x k), co oznacza redukcję
    #       wymiarowości z p do k.
    # --------------------------------------------------------------------------
    transformed_X = X_centered * projection_matrix # Rozmiar (n x k)

    # Obliczenie wyjaśnionej wariancji
    # Proporcja wariancji wyjaśnionej przez i-tą główną składową: λ_i / sum(λ)
    total_variance = sum(eigenvalues)
    explained_variance_ratio = sorted_eigenvalues[1:n_components] ./ total_variance
    
    return transformed_X, explained_variance_ratio
end


"""
    PCAResultAtTime
Struktura do przechowywania wyniku analizy PCA dla pojedynczego punktu w czasie.
"""
struct PCAResultAtTime
    tau::Float64
    transformed_data::Matrix{Float64}
    explained_variance::Vector{Float64}
end


"""
    analyze_pca_over_time(sim_result::modHydroSim.SimResult, sample_times::AbstractVector; n_components=2)
Przeprowadza analizę PCA dla każdej chwili czasu, zwracając wektor wyników.
"""
function analyze_pca_over_time(sim_result::modHydroSim.SimResult, sample_times::AbstractVector; n_components=2)
    solutions = sim_result.solutions
    n_trajectories = length(solutions)
    pca_results_vector = PCAResultAtTime[]

    for tau in sample_times
        data_at_tau = zeros(Float64, n_trajectories, 2)
        for i in 1:n_trajectories
            data_at_tau[i, :] = solutions[i](tau)
        end
        
        transformed_data, explained_variance = pca_math(data_at_tau, n_components)
        result = PCAResultAtTime(tau, transformed_data, explained_variance)
        push!(pca_results_vector, result)
    end
    
    return pca_results_vector
end


"""
    plot_all_pca_results(pca_results::Vector{PCAResultAtTime}, sim_result::modHydroSim.SimResult)
Generuje serię wykresów dla wszystkich dostępnych wyników analizy PCA.
"""
function plot_all_pca_results(pca_results::Vector{PCAResultAtTime}, sim_result::modHydroSim.SimResult)
    initial_temps = [sol.u[1][1] for sol in sim_result.solutions]
    
    for result in pca_results
        total_explained_var = sum(result.explained_variance) * 100
        
        p = scatter(
            result.transformed_data[:, 1], 
            result.transformed_data[:, 2],
            marker_z=initial_temps,
            title="PCA dla τ = $(round(result.tau, digits=2)) (Wariancja: $(round(total_explained_var, digits=2))%)",
            xlabel="Komponent PCA 1",
            ylabel="Komponent PCA 2",
            label="",
            markersize=5,
            alpha=0.8,
            color=:plasma,
            colorbar_title="  Temp. początkowa [MeV]"
        )
        display(p)
    end
end

function plot_pca_at_time(pca_results::Vector{PCAResultAtTime}, sim_result::modHydroSim.SimResult, target_tau::Float64)
    idx = findfirst(res -> res.tau ≈ target_tau, pca_results)
    
    if isnothing(idx)
        available_times = [round(r.tau, digits=3) for r in pca_results]
        error("Nie znaleziono wyników dla tau ≈ $target_tau. Dostępne czasy: $available_times")
    end
    
    result = pca_results[idx]
    initial_temps = [sol.u[1][1] for sol in sim_result.solutions]
    total_explained_var = sum(result.explained_variance) * 100
    
    p = scatter(
        result.transformed_data[:, 1], 
        result.transformed_data[:, 2],
        marker_z=initial_temps,
        title="PCA dla τ = $(round(result.tau, digits=2)) (Wariancja: $(round(total_explained_var, digits=2))%)",
        xlabel="Komponent PCA 1",
        ylabel="Komponent PCA 2",
        label="",
        markersize=5,
        alpha=0.8,
        color=:plasma,
        colorbar_title="  Temp. początkowa [MeV]"
    )
    display(p)
end

end # --- koniec modułu PCAWorkflow ---


#    julia> include("pca.jl")
#

#    julia> using .PCAWorkflow
#

#    julia> settings = SimSettings(n_points=300, tspan=(0.2, 1.2))
#    julia> sim_data = run_simulation(settings)
#

#    julia> czasy_wykres = [0.2, 0.6, 1.0]
#    julia> wyniki_analizy = analyze_pca_over_time(sim_data, czasy_wykres)

#    julia> plot_pca_at_time(wyniki_analizy, sim_data, 0.6)

#    julia> plot_all_pca_results(wyniki_analizy, sim_data)
# --------------------------------------------------------------------------
