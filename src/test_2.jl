using HydroAttractors
using Plots

function kadr(simres, t)
    (Ts, As) = TA(simres.solutions,t)
    p = plot(title="Kadr w przestrzeni fazowej, t = $t",
             xlabel="Temperatura T [MeV]",
             ylabel="Anizotropia A")
    plot!(p, Ts, As, seriestype=:scatter, label="")
    display(p)
end


println("Biblioteki wczytane. Program Do robienia wykresu na bazie run_simulation z lib.jl")

ustawienia = SimSettings(tspan=(0.2, 1.0))

wynik_symulacji = run_simulation(ustawienia)

println("Symulacja zakończona. Generowanie wykresu...")

kadr(wynik_symulacji, 0.5)

println("\nWykres gotowy. Naciśnij Enter, aby zakończyć.")

readline()
