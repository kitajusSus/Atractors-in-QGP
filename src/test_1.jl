include("lib.jl")
using .modHydroSim
using Plots
gr() 
println("Biblioteki wczytane. Program Do robienia wykresu na bazie run_simulation z lib.jl")


ustawienia = SimSettings()


parametry = PARAMS_SYM
wynik_symulacji = run_simulation(parametry, ustawienia)


println("Symulacja zakończona. Generowanie wykresu...")


p = plot(
    title="Trajektorie w przestrzeni fazowej",
    xlabel="Temperatura T [MeV]",
    ylabel="Anizotropia A"
)

# Iterujemy przez KAŻDE rozwiązanie zwrócone przez `run_simulation`
for sol in wynik_symulacji.solutions
    #  `plot!` (z wykrzyknikiem), POZWALA DODAĆ linię do istniejącego wykresu `p`.
    # `idxs=2` wybiera drugą zmienną (A) do narysowania.
    # `label=""` 
    plot!(p, sol, idxs=(1,2), label="")
end


# --- KROK 4: Wyświetlenie i zatrzymanie skryptu ---
display(p)
println("\nWykres gotowy. Naciśnij Enter, aby zakończyć.")
readline()
