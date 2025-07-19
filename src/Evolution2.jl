# 1. Importowanie potrzebnych pakietów
using CSV
using DataFrames
using Plots
using Printf

# Nazwa pliku wejściowego z danymi
const input_filename = "wyniki_symulacji_QGP2.csv"
const output_file = "evolution2_ruch_punktow.gif"

# Sprawdzenie, czy plik z danymi istnieje
if !isfile(input_filename)
    println("BŁĄD: Nie znaleziono pliku '$input_filename'.")
    println("Uruchom najpierw skrypt generujący dane, aby go utworzyć.")
    exit()
end

# 2. Wczytanie danych z pliku CSV
println("Wczytywanie danych z pliku '$input_filename'...")
df = CSV.read(input_filename, DataFrame)
println("Dane wczytane. Znaleziono $(nrow(df)) wierszy.")

# 3. Przygotowanie do animacji
# Znajdź wszystkie unikalne kroki czasowe, które będą naszymi klatkami
unique_taus = sort(unique(df.tau))
n_frames = length(unique_taus)
println("Znaleziono $n_frames klatek do animacji.")

# Parametry animacji i wykresu
const animation_fps = 20
const plot_xlims = (50,400 )
const plot_ylims = (0, 4.9)

# 4. Główna pętla tworząca animację
println("Rozpoczynam generowanie animacji ruchu punktów...")

# Użyjemy ciemnego motywu dla lepszego kontrastu
theme(:dark) 

anim = @animate for (i, τ_current) in enumerate(unique_taus)
    
    # Wyświetlanie postępu w terminalu
    print("\rGenerowanie klatki $i / $n_frames (τ = $(round(τ_current, digits=2)))")

    # Filtrujemy cały DataFrame, aby uzyskać wszystkie punkty dla bieżącego τ
    current_frame_data = filter(row -> row.tau == τ_current, df)
    
    # Tworzymy wykres punktowy dla tej jednej klatki
    scatter(current_frame_data.T_at_tau, current_frame_data.A_at_tau.*1000,
        # Opcje wykresu
        title = @sprintf("Chmura punktów w przestrzeni stanów (τ = %.2f fm/c)", τ_current),
        xlabel = "Temperatura T [MeV]",
        ylabel = "Anizotropia A",
        xlims = plot_xlims,
        ylims = plot_ylims,
        legend = false,
        framestyle = :box,
        background_color = :black,
        
        # Opcje punktów
        markersize = 3,
        markerstrokewidth = 0.2, # Dodajemy lekką obwódkę dla czytelności
        markerstrokecolor = :white,
        alpha = 0.8,
        
        # Opcje kolorowania
        zcolor = current_frame_data.T_0, # Kolor zależy od początkowej temperatury
        c = :plasma, # Paleta kolorów
        colorbar_title = "\n  T₀ [MeV]"
    )
end

# 5. Zapisywanie animacji do pliku GIF
println("Zapisywanie animacji do pliku '$output_file'. Może to chwilę potrwać...")
gif(anim, output_file, fps = animation_fps)

println("Gotowe! Animacja ruchu punktów została pomyślnie utworzona.")
