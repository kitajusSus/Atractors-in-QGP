using AttractorsQGP
using StaticArrays
using Printf

println("Generating matched datasets for MIS and HJSW...")

# Parametry
n_points = 1000
tspan = (0.2, 10.0)
T_range = (100.0, 1000.0)
A_range = (-1.0, 6.0)
seed = 5
saveat = 0.05

model_mis = MISModel()
model_hjsw = HJSWModel()

# 1. Wygeneruj warunki początkowe MIS [T0, A0]
println("1. Generating initial conditions for MIS (n_points=$n_points, seed=$seed)...")
ics_mis = generate_initial_conditions(
    n_points;
    T_range = T_range,
    A_range = A_range,
    temperature_unit = :MeV,
    seed = seed
)

# 2. Symuluj trajektorie MIS
println("2. Simulating MIS trajectories...")
sols_mis = generate_trajectories(model_mis, ics_mis, tspan; saveat = saveat)
dataset_mis = build_dataset(sols_mis)

output_mis = "datasets/hjsw_lpca/mis_matched_seed5.h5"
save_dataset(output_mis, dataset_mis)
println("   Saved MIS dataset to: $output_mis (size: $(size(dataset_mis)))")

# 3. Przygotuj pasujące 3D warunki początkowe dla HJSW [T0, A0, B0]
println("3. Computing matching initial conditions [T0, A0, B0] for HJSW...")
tau0 = tspan[1]
ics_hjsw = Vector{SVector{3, Float64}}(undef, length(ics_mis))

for i in eachindex(ics_mis)
    ic_mis = ics_mis[i] # [T0_fm, A0]
    du_mis = AttractorsQGP.rhs(ic_mis, model_mis, tau0)
    dA0_mis = du_mis[2]
    B0 = tau0 * dA0_mis # B = tau * dA/dtau
    ics_hjsw[i] = SVector{3, Float64}(ic_mis[1], ic_mis[2], B0)
end

# 4. Symuluj trajektorie HJSW
println("4. Simulating HJSW trajectories with matched initial conditions...")
sols_hjsw = generate_trajectories(model_hjsw, ics_hjsw, tspan; saveat = saveat)
dataset_hjsw = build_dataset(sols_hjsw)

output_hjsw = "datasets/hjsw_lpca/hjsw_matched_seed5.h5"
save_dataset(output_hjsw, dataset_hjsw)
println("   Saved HJSW dataset to: $output_hjsw (size: $(size(dataset_hjsw)))")

println("Done generating matched datasets!")
