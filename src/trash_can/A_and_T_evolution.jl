# ===================================================================
#      GENERATOR ANIMACJI Z DANYCH SYMULACYJNYCH QGP Z PLIKU CSV
# ===================================================================

#w
using CSV
using DataFrames
using Plots
using Printf

# Nazwa pliku wejściowego z danymi
const input_filename = "dane_log_tau.csv"
const output_file = "animacja_z_danych.gif"

# Sprawdzenie, czy plik z danymi istnieje
if !isfile(input_filename)
    println("BŁĄD: Nie znaleziono pliku '$input_filename'.")
    println("Upewnij się, że plik z danymi znajduje się w tym samym folderze co skrypt.")
    println("Możesz go wygenerować za pomocą poprzedniego programu.")
    exit()
end

# 2. Wczytanie danych z pliku CSV
println("Wczytywanie danych z pliku '$input_filename'...")
df = CSV.read(input_filename, DataFrame)
println("Dane wczytane pomyślnie. Znaleziono $(nrow(df)) wierszy.")

# 3. Przygotowanie do animacji
# Znajdź wszystkie unikalne kroki czasowe, które będą naszymi klatkami
unique_taus = sort(unique(df.tau))
n_frames = length(unique_taus)
println("Znaleziono $n_frames unikalnych kroków czasowych do animacji.")

# Grupujemy dane po identyfikatorze symulacji (Run_ID)
# To pozwoli nam łatwo rysować poszczególne trajektorie.
grouped_data = groupby(df, :Run_ID)

# Parametry animacji i wykresu
const animation_fps = 20
const plot_xlims = (150, 800)
const plot_ylims = (-1.5, 5.5)

# 4. Główna pętla tworząca animację
println("Rozpoczynam generowanie animacji...")
theme(:dark) # Ciemny motyw dla lepszego kontrastu

anim = @animate for (i, τ_current) in enumerate(unique_taus)
    
    print("\rGenerowanie klatki $i / $n_frames (τ = $(round(τ_current, digits=2)))")

    # Inicjalizujemy pusty wykres dla każdej klatki
    p = plot(
        title = @sprintf("Ewolucja do atraktora (τ = %.2f fm/c)", τ_current),
        xlabel = "Temperatura T [MeV]",
        ylabel = "Anizotropia A",
        xlims = plot_xlims,
        ylims = plot_ylims,
        legend = false,
        framestyle = :box,
        #background_color = :black,
        #fg_color_legend= :white,
        #fg_color_text=:white,
        #fg_color_axis=:white
    )

    # --- Rysowanie trajektorii ("ogonów") ---
    # Iterujemy po każdej grupie (czyli po każdej pojedynczej symulacji)
    for group in grouped_data
        # Filtrujemy dane dla bieżącej trajektorii do aktualnego momentu czasu
        path_data = filter(row -> row.tau <= τ_current, group)
        
        if !isempty(path_data)
            # Kolor jest stały dla całej trajektorii i zależy od T_0
            line_color_value = path_data.T_0[1] 
            
            plot!(p, path_data.T_at_tau, path_data.A_at_tau,
                linewidth=1,
                alpha=1, # "Ogon" jest lekko przezroczysty
                line_z = line_color_value,
                c = :plasma,
            )
        end
    end
    
    # --- Rysowanie aktualnych pozycji ("główek") ---
    # Filtrujemy cały DataFrame, aby uzyskać wszystkie punkty dla bieżącego τ
    head_data = filter(row -> row.tau == τ_current, df)
    
    scatter!(p, head_data.T_at_tau, head_data.A_at_tau,
        markersize = 2.5,
        markerstrokewidth = 0,
        zcolor = head_data.T_0, # Kolor punktu na podstawie T_0
        c = :plasma,
        colorbar_title = "\n  T₀ [MeV]"
    )
    
    p # Zwracamy gotowy wykres jako klatkę animacji
end

# 5. Zapisywanie animacji do pliku GI:F
println("\n\nZapisywanie animacji do pliku '$output_file'. Może to chwilę potrwać...")
gif(anim, output_file, fps = animation_fps)

println("Gotowe! Animacja została pomyślnie utworzona.")
