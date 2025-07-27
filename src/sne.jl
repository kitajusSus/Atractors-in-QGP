include("lib.jl")
using .modHydroSim
using Plots
"""
Cała implementacja bazuje na pracy Hinton'a i Roweis'a z 2003 roku: 
URL : https://www.cs.toronto.edu/~hinton/absps/sne.pdf
*Stochastic Neighbor Embedding*

"""

module modSNE

using Statistics
using Random
using LinearAlgebra
using Printf

export SNESettings, run_sne

"""
    SNESettings(output_dims, perplexity, iterations, ...)

Struktura przechowująca wszystkie hiperparametry dla algorytmu SNE.
"""
struct SNESettings
    output_dims::Int
    perplexity::Float64
    iterations::Int
    learning_rate::Float64
    initial_jitter::Float64
end

function SNESettings(; 
                     output_dims=2, 
                     perplexity=20.0, 
                     iterations=1000, 
                     learning_rate=100.0, 
                     initial_jitter=0.5
                     )  
  # Wywołuje domyślny, wewnętrzny konstruktor
    return SNESettings(output_dims, perplexity, iterations, learning_rate, initial_jitter)
end

"""
Oblicza macierz kwadratów odległości Euklidesowych między wszystkimi parami punktów.
"""
function calculate_squared_distances(X::Matrix{Float64})
    (n_points, _) = size(X)
    dist_sq = zeros(Float64, n_points, n_points)
    for i in 1:n_points
        for j in (i+1):n_points
            d = norm(X[i, :] - X[j, :])
            dist_sq[i, j] = dist_sq[j, i] = d^2
        end
    end
    return dist_sq
end


"""
Dla danego `σ` (sigma) i wektora odległości, oblicza rozkład prawdopodobieństwa
`P` oraz jego entropię.
"""
function calculate_p_and_entropy(dist_row_sq::Vector{Float64}, sigma::Float64, point_idx::Int)
    # Wzór (1) i (2) z pracy -sne hinton -   exp(-||xi-xj||^2 / 2σi_i^2)
    p = exp.(-dist_row_sq / (2 * sigma^2))
    p[point_idx] = 0.0 # Prawdopodobieństwo do samego siebie jest zerowe
    
    sum_p = sum(p)
    if sum_p == 0.0 || isnan(sum_p)
        return (fill(1e-12, length(p)), 0.0) # dzielenie przez zero
    end

    p_normalized = p ./ sum_p
    
    # Obliczamy entropię Shannona H = -Σ p log2(p)
    # uzywamy logarytmu naturalnego, wiec Perplexity = exp(H)
    entropy = -sum(p_normalized .* log.(max.(p_normalized, 1e-12))) # max dla stabilnosci num.

    return p_normalized, entropy
end

"""
Przeprowadza wyszukiwanie binarne, aby znaleźć `σ` (sigma), które
daje zadaną perpleksję dla jednego punktu.
"""
function find_sigma_for_perplexity(dist_row_sq::Vector{Float64}, perplexity::Float64, point_idx::Int)
    target_entropy = log(perplexity)
    
    sigma_min = 0.0
    sigma_max = Inf
    sigma_val = 1.0 # Początkowy strzał

    # Wyszukiwanie binarne przez 50 kroków
    for _ in 1:50
        (_, current_entropy) = calculate_p_and_entropy(dist_row_sq, sigma_val, point_idx)
        
        if current_entropy < target_entropy
            sigma_min = sigma_val
            sigma_val = (isinf(sigma_max)) ? sigma_val * 2 : (sigma_val + sigma_max) / 2.0
        else
            sigma_max = sigma_val
            sigma_val = (sigma_val + sigma_min) / 2.0
        end
    end
    return sigma_val
end

"""
Oblicza macierz `P` prawdopodobieństw wysokowymiarowych `p_j|i`.
Jest to najbardziej kosztowna obliczeniowo część inicjalizacji.
"""
function calculate_P_matrix(X::Matrix{Float64}, perplexity::Float64)
    (n_points, _) = size(X)
    dist_sq = calculate_squared_distances(X)
    P = zeros(Float64, n_points, n_points)

    println("Obliczanie macierzy P (prawdopodobieństw wysokowymiarowych)...")
    for i in 1:n_points
        print("\rPrzetwarzanie punktu $i/$n_points...")
        # Znajdź optymalne sigma dla punktu i
        sigma_i = find_sigma_for_perplexity(dist_sq[i, :], perplexity, i)
        # Oblicz wiersz macierzy P dla tego sigma
        (p_row, _) = calculate_p_and_entropy(dist_sq[i, :], sigma_i, i)
        P[i, :] = p_row
    end
    println("\nObliczanie macierzy P zakończone.")
    return P
end


"""
Oblicza macierz `Q` prawdopodobieństw niskowymiarowych `q_j|i`.
Wymaga ponownego obliczania w każdej iteracji.
Zgodnie z pracą, wariancja w przestrzeni niskowymiarowej jest stała (1/2).
"""
function calculate_Q_matrix(Y::Matrix{Float64})
    (n_points, _) = size(Y)
    dist_sq = calculate_squared_distances(Y)
    
    # Wzór (3) z pracy: q_j|i = exp(-||yi-yj||^2) / Σk exp(-||yi-yk||^2)
    # wariancja jest 1/2, więc 2σ^2 = 1
    q_unnormalized = exp.(-dist_sq)
    # Ustawiamy przekątną na zero
    q_unnormalized[diagind(q_unnormalized)] .= 0.0
    
    row_sums = sum(q_unnormalized, dims=2)
    Q = q_unnormalized ./ max.(row_sums, 1e-12) # max dla stabilnosci num.

    return Q
end


"""
Oblicza gradient funkcji kosztu `C` względem położeń `y_i`.
Implementuje wzór (5) z pracy.
"""
function calculate_gradient(P::Matrix{Float64}, Q::Matrix{Float64}, Y::Matrix{Float64})
    (n_points, dims) = size(Y)
    grad = zeros(Float64, n_points, dims)

    for i in 1:n_points
        for j in 1:n_points
            if i == j continue end
            # Wzór (5): ∂C/∂yi = Σj 2 * (yi - yj) * (p_j|i - q_j|i + p_i|j - q_i|j)
            force_magnitude = 2 * ( (P[j, i] - Q[j, i]) + (P[i, j] - Q[i, j]) )
            grad[i, :] += force_magnitude * (Y[i, :] - Y[j, :])
        end
    end
    return grad
end

"""
Oblicza koszt (sumę dywergencji KL), zgodnie ze wzorem (4).
"""
function calculate_cost(P::Matrix{Float64}, Q::Matrix{Float64})
    # Wzór (4): C = Σi Σj p_j|i log(p_j|i / q_j|i)
    # Dodajemy 1e-12 dla stabilności numerycznej
    cost_matrix = P .* log.((P .+ 1e-12) ./ (Q .+ 1e-12))
    return sum(cost_matrix)
end


"""
    run_sne(X::Matrix{Float64}, settings::SNESettings) -> Matrix{Float64}

Główna funkcja uruchamiająca algorytm SNE.
"""
function run_sne(X::Matrix{Float64}, settings::SNESettings)
    (n_points, _) = size(X)
    
    # Oblicz macierz P (tylko raz)
    P = calculate_P_matrix(X, settings.perplexity)
    
    # Krok 1: Inicjalizacja położeń Y w przestrzeni niskowymiarowej
    # Zgodnie z pracą: losowe lokalizacje blisko zera
    rng = MersenneTwister(123)
    Y = randn(rng, n_points, settings.output_dims) * 0.0001
    
    println("Rozpoczynanie optymalizacji SNE...")
    
    # Główna pętla optymalizacji
    for iter in 1:settings.iterations
        # Krok 2: Oblicz macierz Q dla aktualnych położeń Y
        Q = calculate_Q_matrix(Y)
        
        # Krok 3: Oblicz gradient
        grad = calculate_gradient(P, Q, Y)
        
        # Krok 4: Zaktualizuj położenia Y (spadek gradientu)
        Y .-= settings.learning_rate .* grad
        
        # Krok 5: Dodaj wygaszany szum (annealed noise / jitter)
        # Zgodnie z pracą, szum jest redukowany w czasie
        current_jitter = settings.initial_jitter * (1.0 - iter / settings.iterations)
        if current_jitter > 0
            Y .+= randn(rng, size(Y)) * current_jitter
        end

        # Krok 6: Wyśrodkuje dane
        Y .-= mean(Y, dims=1)
        
        # Raportowanie postępów
        if iter % 50 == 0
            cost = calculate_cost(P, Q)
            @printf("Iteracja %4d/%d: Koszt = %.4f\n", iter, settings.iterations, cost)
        end
    end
    
    println("Optymalizacja SNE zakończona.")
    return Y
end


end # koniec modułu modSNE

# ##############################################################################
# SEKCJA 2: URUCHOMIENIE EKSPERYMENTU
# ##############################################################################

"""
    prepare_hydro_data(sim_settings::SimSettings, n_samples::Int) -> Matrix, Vector

Uruchamia symulację hydrodynamiki i przekształca trajektorie w dane
wysokowymiarowe gotowe do użycia przez SNE.
"""
function prepare_hydro_data(sim_settings::SimSettings, n_time_samples::Int)
    println("--- Krok 1: Generowanie danych z symulacji hydrodynamiki... ---")
    wynik_symulacji = run_simulation(sim_settings)
    n_trajectories = length(wynik_symulacji.solutions)
    
    # Definiujemy czasy, w których będziemy próbkować każdą trajektorię
    t_start, t_end = wynik_symulacji.settings.tspan
    sample_times = range(t_start, stop=t_end, length=n_time_samples)
    
    # Każda trajektoria stanie się jednym wierszem w macierzy X
    # Wymiar D = n_time_samples * 2 (dla T i A)
    D = n_time_samples * 2
    X = zeros(Float64, n_trajectories, D)
    
    # Zapisujemy też początkowe temperatury do późniejszego kolorowania
    initial_temperatures = zeros(Float64, n_trajectories)

    println("--- Krok 2: Przekształcanie trajektorii w wektory wysokowymiarowe... ---")
    for i in 1:n_trajectories
        sol = wynik_symulacji.solutions[i]
        sampled_points = sol(sample_times)
        
        # "Spłaszczamy" próbki [T1, A1, T2, A2, ...] do jednego wektora
        flat_vector = vcat([[p[1], p[2]] for p in sampled_points]...)
        X[i, :] = flat_vector
        
    initial_temperatures[i] = sol.u[1][1] # Temperatura początkowa
    end
    
    println("Przygotowanie danych zakończone. Wymiar danych wejściowych: $n_trajectories x $D")
    return X, initial_temperatures
end



 
#Dodatkowa funkcja 

function main()
    hydro_settings = SimSettings(
        n_points=250, 
        tspan=(0.2, 0.8)
    )
    n_time_samples = 5 # Liczba próbek z każdej trajektorii

    # --- Ustawienia algorytmu SNE ---
    sne_settings = modSNE.SNESettings(
        output_dims=2,
        perplexity=25,
        iterations=2,
        learning_rate=100,
        initial_jitter=0.5
    )

    # Krok 1 & 2: Przygotuj dane
    X, initial_temps = prepare_hydro_data(hydro_settings, n_time_samples)

    # Krok 3: Uruchom SNE
    println("\n--- Krok 3: Uruchamianie algorytmu SNE... ---")
    Y = modSNE.run_sne(X, sne_settings)

    # Krok 4: Wizualizacja
    println("\n--- Krok 4: Generowanie wykresu wynikowego... ---")
    p = scatter(
        Y[:, 1], Y[:, 2],
        marker_z=initial_temps, # Koloruj punkty wg. temperatury początkowej
        title="Wizualizacja SNE trajektorii hydrodynamiki",
        xlabel="Komponent SNE 1",
        ylabel="Komponent SNE 2",
        label="",
        markersize=5,
        alpha=0.8,
        color=:plasma,
        colorbar_title="  Temp. początkowa [MeV]"
    )
    
    display(p)
    println("Gotowe. Wyświetlono wykres.")
end

# Uruchomienie głównej funkcji
main()

