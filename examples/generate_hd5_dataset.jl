using DifferentialEquations
using OrdinaryDiffEqRosenbrock
using HDF5
using Random
using StaticArrays

# Stała przeliczeniowa: T [fm^-1] = T [MeV] / HBARC_MEV_FM
const HBARC_MEV_FM = 197.3269804

# Parametry modelu HJSW (N=4 SYM)
const Ω_R = 9.8
const Ω_I = 8.629
const C_η = 1.0 / (4.0 * π)
const Ω2 = Ω_R^2 + Ω_I^2

"""
Równania różniczkowe ewolucji modelu HJSW w przepływie Bjorkena:
u = [T, A, B]
t = czas własny tau [fm/c]
"""
function hjsw_rhs(u, p, t)
    T, A, B = u[1], u[2], u[3]

    dT = (1.0 / 18.0) * (A - 6.0) * T / t
    dA = B / t

    term1 = A * (-(11.0 * B) / 18.0 - t^2 * T^2 * Ω2)
    term2 = -(4.0 / 27.0) * (A^2) * (3.0 * t * Ω_I * T - 1.0)
    term3 = -(1.0 / 27.0) * (A^3)
    term4 = B * ((2.0 / 3.0) - 2.0 * t * Ω_I * T)
    term5 = 8.0 * C_η * t * T * Ω2

    dB = (term1 + term2 + term3 + term4 + term5) / t

    return SVector(dT, dA, dB)
end

function generate_hjsw_hd5_dataset(;
        n_runs::Int = 1000,
        output_file::String = "datasets/hjsw_dataset.hd5",
        tspan::Tuple{Float64, Float64} = (0.2, 3.0),
        dt_save::Float64 = 0.01,
        T_range_MeV::Tuple{Float64, Float64} = (100.0, 1000.0),
        A_range::Tuple{Float64, Float64} = (-1.0, 6.0),
        B_range::Tuple{Float64, Float64} = (-1.0, 1.0),
        seed::Int = 42
    )
    rng = Xoshiro(seed)

    # Zakres T w jednostkach 1/fm (fermi=1)
    T_min_fm = T_range_MeV[1] / HBARC_MEV_FM
    T_max_fm = T_range_MeV[2] / HBARC_MEV_FM

    mkpath(dirname(output_file))

    println("Generowanie $(n_runs) rozwiązań równań HJSW...")
    println("Zakres początkowy T: $(T_range_MeV[1]) MeV ... $(T_range_MeV[2]) MeV ($(round(T_min_fm, digits=4)) ... $(round(T_max_fm, digits=4)) 1/fm)")
    println("Zakres początkowy A: $(A_range[1]) ... $(A_range[2])")
    println("Zakres początkowy B: $(B_range[1]) ... $(B_range[2])")
    println("Zakres czasu tau: $(tspan[1]) ... $(tspan[2]) fm/c (krok = $(dt_save))")

    h5open(output_file, "w") do h5file
        saved_count = 0
        for i in 1:n_runs
            # Losowanie warunków początkowych
            T0 = rand(rng) * (T_max_fm - T_min_fm) + T_min_fm
            A0 = rand(rng) * (A_range[2] - A_range[1]) + A_range[1]
            B0 = rand(rng) * (B_range[2] - B_range[1]) + B_range[1]

            u0 = SVector{3, Float64}(T0, A0, B0)

            prob = ODEProblem(hjsw_rhs, u0, tspan)
            sol = solve(prob, Rodas5(), saveat = dt_save, reltol = 1e-8, abstol = 1e-8)

            if sol.retcode == ReturnCode.Success
                # Tworzenie macierzy [t, T, A, B]
                n_pts = length(sol.t)
                matrix = Matrix{Float64}(undef, n_pts, 4)
                for k in 1:n_pts
                    matrix[k, 1] = sol.t[k]
                    matrix[k, 2] = sol.u[k][1]
                    matrix[k, 3] = sol.u[k][2]
                    matrix[k, 4] = sol.u[k][3]
                end

                # Zapis do zestawu danych w HDF5 (np. "run-1", "run-2", ...)
                h5file["run-$i"] = matrix
                saved_count += 1
            end

            if i % 200 == 0 || i == n_runs
                println("Przetworzono $i / $n_runs rozwiązań...")
            end
        end
        println("Zapisano $saved_count rozwiązań do pliku: $output_file")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    generate_hjsw_hd5_dataset()
end
