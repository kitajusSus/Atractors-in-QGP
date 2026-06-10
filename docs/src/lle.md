By dostać wykresy analizy spectrum widma wartości własnych robiłem w repl takie coś 


- `mis = load_dataset("datasets/....")`
- `taus` => `[0.2, 0.25, 0.35, 0.45,0.55,0.65]`



```julia
for (i, tau) in enumerate(taus)
           fig = plot_lle_spectrum_statistics(mis; τ = tau, k_values = 50:5:250, ile_λ = 3)
           filename = "plots/lle/lle_eigen2_tau_$(tau).png"
           save(filename, fig)
       end

for (i, tau) in enumerate(taus)
           fig = plot_lle_spectrum_statistics(mis; τ = tau, k_values = 5:5:50, ile_λ = 3)
           filename = "plots/lle/lle_eigen_tau_$(tau).png"
           save(filename, fig)
       end

```
