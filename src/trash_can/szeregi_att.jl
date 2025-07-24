# -----------------------------------------------------------------
# KROK 1: ZAŁADOWANIE POTRZEBNYCH PAKIETÓW
# -----------------------------------------------------------------
using DifferentialEquations
using Plots
using LaTeXStrings # Pakiet do obsługi notacji LaTeX w etykietach

# -----------------------------------------------------------------
# KROK 2: DEFINICJA PARAMETRÓW FIZYCZNYCH
# -----------------------------------------------------------------
# Używamy wartości dla teorii N=4 SYM (standard w badaniach modelowych)
const C_eta = 1 / (4π)
const C_tau_pi = (2 - log(2)) / (2π)
const C_lambda_1 = 1 / (2π)

const p = (C_eta, C_tau_pi, C_lambda_1)

# -----------------------------------------------------------------
# KROK 3: DEFINICJA RÓWNANIA RÓŻNICZKOWEGO
# -----------------------------------------------------------------
function brsss_attractor!(dA, A_vec, p, w)
    A = A_vec[1]
    C_eta, C_tau_pi, C_lambda_1 = p

    # Zabezpieczenie przed dzieleniem przez zero, jeśli mianownik się wyzeruje
    if abs(1 + (1/12)*A) < 1e-9
        dA[1] = 0.0
        return
    end
    
    licznik = ( (1/3)*C_tau_pi + (C_lambda_1 / (8*C_eta))*w ) * A^2 + (3/2)*w*A - 12*C_eta
    mianownik = C_tau_pi * w * (1 + (1/12)*A)
    
    dA[1] = -licznik / mianownik
end

# -----------------------------------------------------------------
# KROK 4: USTAWIENIA SYMULACJI I TWORZENIE WYKRESU
# -----------------------------------------------------------------
w_span = (0.1, 10.0) # Zakres czasu 'w' dla obliczeń
initial_conditions = -1.0:0.5:5.0 # Zestaw różnych warunków początkowych dla A

# Inicjalizacja pustego wykresu z estetycznymi etykietami i tytułem
p_plot = plot(
    xlabel=L" czas zredukowany, $w = \tau T$",
    ylabel=L"Anizotropia, $\mathcal{A} = (P_T - P_L)/P_{eq}$",
    title="Atraktor Hydrodynamiczny i Rozwinięcia Gradientowe",
    xlims=(0, 10), # Zakres osi X od 0 do 10
    ylims=(-1.5, 5.5),
    legend=:topright,
    framestyle=:box,
    tickfont=font(12),  # Ustawienie rozmiaru czcionki na osiach
    guidefont=font(14), # Ustawienie rozmiaru etykiet osi
    titlefont=font(16)  # Ustawienie rozmiaru tytułu
)

# Pętla rozwiązująca równanie dla każdego warunku początkowego
println("Obliczanie trajektorii dla różnych warunków początkowych...")
for (i, A0) in enumerate(initial_conditions)
    prob = ODEProblem(brsss_attractor!, [A0], w_span, p)
    # Zwiększamy gęstość punktów, aby krzywe były gładsze
    sol = solve(prob, Tsit5(), saveat=0.01) 
    # Dodajemy etykietę tylko dla pierwszej krzywej, aby uniknąć bałaganu w legendzie
    plot!(p_plot, sol, vars=(0,1), color=:dodgerblue, lw=1.5, label= i==1 ? "Trajektorie numeryczne" : "")
end
println("Gotowe.")

# -----------------------------------------------------------------
# KROK 5: OBLICZENIE I NARYSOWANIE SAMEGO ATRAKTORA
# -----------------------------------------------------------------
println("Obliczanie rozwiązania atraktorowego...")
# Zaczynamy bardzo blisko w=0, aby precyzyjnie znaleźć atraktor
w_span_attractor = (0.01, 10.0) 
# Warunek początkowy na atraktorze dla małych w
A_attractor_start = 6 * sqrt(C_eta / C_tau_pi) 

prob_attractor = ODEProblem(brsss_attractor!, [A_attractor_start], w_span_attractor, p)
sol_attractor = solve(prob_attractor, Tsit5(), saveat=0.01)

# Rysujemy atraktor grubą, wyróżniającą się linią
plot!(p_plot, sol_attractor, vars=(0,1), color=:red, lw=3, label=L"Atraktor $\mathcal{A}_{att}$")
println("Gotowe.")

# -----------------------------------------------------------------
# KROK 6: DODANIE ROZWINIĘĆ ASYMPTOTYCZNYCH
# -----------------------------------------------------------------
# Definiujemy zakres 'w', dla którego rysujemy przybliżenia.
# Zaczynamy od w>0, gdzie rozwinięcia mają sens.
w_values = 0.5:0.01:10.0 

# Definicje szeregów asymptotycznych dla BRSSS i MIS do drugiego rzędu
A_brsss(w) = 8*C_eta/w + (16*C_eta*(C_tau_pi - C_lambda_1)) / (3*w^2)
A_mis(w) = 8*C_eta/w + (16*C_eta*C_tau_pi) / (3*w^2)

# Rysujemy obie funkcje na wykresie
plot!(p_plot, w_values, w -> A_brsss(w), # Użycie anonimowej funkcji dla zwięzłości
    color=:magenta,
    linestyle=:dash,
    lw=2.5,
    label="BRSSS (2-gi rząd)"
)

plot!(p_plot, w_values, w -> A_mis(w),
    color=:darkgreen,
    linestyle=:dot,
    lw=2.5,
    label="MIS (2-gi rząd)"
)

# -----------------------------------------------------------------
# KROK 7: ZAPISANIE WYNIKU DO PLIKU
# -----------------------------------------------------------------
nazwa_pliku = "atraktor_hydrodynamiczny_szeregi.png"
# POPRAWIONA LINIA: Dodano brakujący nawias zamykający
savefig(p_plot, nazwa_pliku) 

println("Wykres został zapisany do pliku: ", nazwa_pliku)
