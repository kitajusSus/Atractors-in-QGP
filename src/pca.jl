
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

export SimSettings, PCAResultAtTime, full_analysis, 
       p_pca, p_explained_variance,
       calc_pca

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
    # Krok 5: Transformacja danych do nowej przestrzeni#
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


# --- SEKCJA 2: STRUKTURY I FUNKCJE ANALIZY ---

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
    analyze_pca(sim_result::modHydroSim.SimResult; n_steps=50, n_components=2)

Przeprowadza analizę PCA w `n_steps` punktach czasowych w całym zakresie symulacji.
"""
function calc_pca(sim_result::modHydroSim.SimResult; n_steps=50, n_components=2)
    solutions = sim_result.solutions
    t_start, t_end = sim_result.settings.tspan
    sample_times = range(t_start, t_end, length=n_steps)
    
    n_trajectories = length(solutions)
    pca_results_vector = PCAResultAtTime[]

    for tau in sample_times
        # Zbieranie danych (T, A) dla wszystkich trajektorii w chwili tau
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


# --- SEKCJA 3: WIZUALIZACJA I GŁÓWNY WORKFLOW ---

"""
    p_explained_variance(pca_results::Vector{PCAResultAtTime})

Tworzy wykres pokazujący, jak wariancja wyjaśniona przez PCA
zmienia się w czasie.  wykres na bazie PRL 125
"""
function p_explained_variance(settings::modHydroSim.SimSettings)
    wyniki_symulacji = modHydroSim.run_simulation(settings)
    pca_results = calc_pca(wyniki_symulacji, n_steps=50)
    times = [res.tau for res in pca_results]
    
    # Zakładamy 2 komponenty
    var_pc1 = [res.explained_variance[1] for res in pca_results]
    var_pc2 = [res.explained_variance[2] for res in pca_results]
    println("Wykres dla SimSettings: $(settings)")
    p = plot(
        times, 
        [var_pc1, var_pc2],
        title="Wariancja wyjaśniona przez komponenty PCA w czasie",
        xlabel="Czas τ [fm/c]",
        ylabel="Proporcja wyjaśnionej wariancji",
        label=["Komponent 1" "Komponent 2"],
        legend=:right,
        linewidth=2,
        xticks=[0.2,0.27,0.4,0.61],
    )
    
    hline!(p, [1.0], linestyle=:dash, color=:black, label="100%", alpha=0.5)
    
    return p
end


"""
    p_pca(pca_results::Vector{PCAResultAtTime}, sim_result::modHydroSim.SimResult, target_tau::Float64)

Wyszukuje i wyświetla wykres PCA dla konkretnej chwili czasu `target_tau`.
"""
function p_pca(pca_results::Vector{PCAResultAtTime}, sim_result::modHydroSim.SimResult, target_tau::Float64)
    # Znajdź najbliższy dostępny czas w wynikach
    idx = findmin(res -> abs(res.tau - target_tau), pca_results)[2]
    
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
        xlims=(-1.5, 1.5),
        ylims=(-0.55, 0.51),
        markersize=5,
        alpha=0.8,
        color=:plasma,
        colorbar_title="  Temp. początkowa [MeV]"
    )
    display(p)
    return p
end

"""
    full_analysis(; settings=SimSettings(), n_pca_steps=50)

Uruchamia cały workflow: symulację, analizę PCA, wyświetla kluczowe wyniki
i zwraca artefakty do dalszej analizy.
"""
function full_analysis(; settings=SimSettings(), n_pca_steps=50)
    println("--- Krok 1: Uruchamianie symulacji hydro... ---")
    sim_result = modHydroSim.run_simulation(settings)
    println("Symulacja zakończona. Liczba trajektorii: $(length(sim_result.solutions)).")
    
    println("\n--- Krok 2: Przeprowadzanie analizy PCA w czasie... ---")
    pca_results = calc_pca(sim_result, n_steps=n_pca_steps)
    println("Analiza PCA zakończona dla $n_pca_steps punktów w czasie.")
    
    println("\n--- Krok 3: Generowanie wykresu wyjaśnionej wariancji... ---")
    variance_plot = p_explained_variance(pca_results)
    display(variance_plot)
    
    
    println("Wykres dla koncowej wariancji wyjaśnionej przez PCA gotowy.")
    p_pca(pca_results, sim_result, 0.4)

    final_var_pc1 = pca_results[end].explained_variance[1] * 100
    println("\nAnaliza zakończona. Wariancja wyjaśniona przez PC1 na końcu ewolucji: $(round(final_var_pc1, digits=2))%")
    
    return sim_result, pca_results
end


end # koniec modułu PCAWorkfloexport SimSettings, run_simulation, analyzePCA, p_allPCA, p_PCAex, pca_mathw
