using GLMakie
using LaTeXStrings
using Statistics
using LinearAlgebra
using AtractorsQGP

function normalizuj_max(X_features)
    T_col = X_features[:, 1]
    A_col = X_features[:, 2]
    
    T_max = maximum(abs.(T_col))
    A_max = maximum(abs.(A_col))
    
    T_norm = T_col ./ (T_max == 0 ? 1.0 : T_max)
    A_norm = A_col ./ (A_max == 0 ? 1.0 : A_max)
    
    return hcat(T_norm, A_norm)
end

function get_tau_slice(dataset::AbstractMatrix{<:Real}, tau::Real; atol::Real=1e-12, feature_cols::AbstractVector{<:Integer}=collect(2:size(dataset, 2)))
    idx = findall(isapprox.(view(dataset, :, 1), tau; atol=atol))
    Xtau = Matrix{Float64}(dataset[idx, feature_cols])
    return idx, Xtau
end

include("lpca.jl")

function ex_lpca(dataset_loadet, k_bazowe::Int, n_slices::Int)
    taus = sort(unique(dataset_loadet[:, 1]))
    step = max(1, length(taus) ÷ n_slices)
    wybrane_tau = taus[1:step:end][1:n_slices]

    k_zmienione = round(Int, k_bazowe * 2)
    tolerancja = 0.01

    println("Tau | Śr. wymiar dla K=$(k_bazowe) | Śr. wymiar dla K=$(k_zmienione) | Zmiana (%)")
    println("-" ^ 75)
    procenty = Float64[]
    tau_values = Float64[]
    mean_k1_values = Float64[]
    mean_k2_values = Float64[]

    for tau in wybrane_tau
        idx, X_tau = get_tau_slice(dataset_loadet, tau; feature_cols=[2, 3])
        X_norm = normalizuj_max(X_tau)
        
        dims_k1 = dims(X_norm; k=k_bazowe, tol=tolerancja)
        dims_k2 = dims(X_norm; k=k_zmienione, tol=tolerancja)
        
        mean_k1 = mean(dims_k1)
        mean_k2 = mean(dims_k2)
        
        diff_percent = mean_k1 > 0 ? ((mean_k2 - mean_k1) / mean_k1) * 100 : 0.0
        

        println(tau, " | ", mean_k1, " | ", mean_k2, " | ", diff_percent, "%")
        push!(procenty, diff_percent)
        push!(tau_values, tau)
        push!(mean_k1_values, mean_k1)
        push!(mean_k2_values, mean_k2)
    end
    return tau_values, mean_k1_values, mean_k2_values, procenty;
end


function main_local_pca(
    dataset_loaded::AbstractMatrix{<:Real};
    n_slices::Int = 15,
    tablica_k::Vector{Int} = [10, 20, 30, 40]
)
    set_publication_theme()

    fig = Figure(size = (950, 600))
    ax = Axis(
        fig[1, 1],
        title = L"\text{Zależność wymiaru lokalnego od czasu własnego } \tau",
        xlabel = L"\tau\,[\mathrm{fm}/c]",
        ylabel = L"\text{Średni wymiar lokalny}",
        limits = (nothing, nothing, 0.8, 2.2)
    )

    palette = Makie.wong_colors()

    for (i, k_bazowe) in enumerate(tablica_k)
        k_zmienione = k_bazowe * 2
        
        tau_vals, mean_k1, mean_k2, _ = ex_lpca(dataset_loaded, k_bazowe, n_slices)
        
        c = palette[mod1(i, length(palette))]

        band!(
            ax, 
            tau_vals, 
            mean_k1, 
            mean_k2; 
            color = (c, 0.25)
        )

        lines!(
            ax, 
            tau_vals, 
            mean_k1; 
            linewidth = 2.5, 
            color = c, 
            label = L"K \in [%$(k_bazowe), %$(k_zmienione)]"
        )
        
        lines!(
            ax, 
            tau_vals, 
            mean_k2; 
            linewidth = 1.0, 
            color = c, 
            linestyle = :dash
        )
    end

    axislegend(ax, position = :rt)
    
    return fig
end
