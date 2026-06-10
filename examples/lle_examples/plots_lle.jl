function plot_examples_lle(X, Y, labels)
    fig = Figure(size = (1200, 600))
    sizeX = size(X)
    ax_3d = Axis3(fig[1, 1], title = "Oryginalna rozmaitość X $sizeX", azimuth = 0.22 * π)
    scatter!(ax_3d, X[1, :], X[2, :], X[3, :], color = labels, colormap = :jet, markersize = 8)
    # println(Y)
    ax_2d = Axis3(fig[1, 2], title = "Zredukowana przestrzeń Y $(size(Y))", azimuth = 0.22 * π)
    scatter!(ax_2d, Y[1, :], Y[2, :], Y[3, :], color = labels, colormap = :jet, markersize = 8)

    return fig
end
function plot_lle_spectrum_heatmap(dataset; τs, k = 20)
    set_publication_theme()

    spectra_matrix = []

    for τ in τs
        _, Xτ = get_tau_slice(dataset, τ)
        X = Matrix{Float64}(Xτ)'

        if size(X, 2) > k + 2
            λ = lle_spectrum(X; k = k)
            push!(spectra_matrix, λ)
        end
    end

    min_len = minimum(length.(spectra_matrix))
    S = hcat([s[1:min_len] for s in spectra_matrix]...)

    fig = Figure(size = (900, 600))
    ax = Axis(
        fig[1, 1],
        title = L"\text{LLE spectrum evolution over } \tau",
        xlabel = L"\text{eigen index}",
        ylabel = L"\tau"
    )

    heatmap!(ax, S)

    return fig
end

function plot_lle_spectrum_stability(dataset; τ, k_values = 5:5:50)
    set_publication_theme()

    spectra, ks = compute_lle_spectra_over_k(dataset, τ; k_values = k_values)

    # stability = variance of first eigenvalue (or gap)
    stability = Float64[]

    for s in spectra
        push!(stability, s[2] - s[1])  # spectral gap
    end

    fig = Figure(size = (900, 500))
    ax = Axis(
        fig[1, 1],
        title = L"\text{LLE spectral gap stability vs k, } \tau=%$τ",
        xlabel = L"k",
        ylabel = L"\lambda_2 - \lambda_1"
    )

    lines!(ax, ks, stability, linewidth = 3)

    return fig
end

function plot_examples_lle_3d(X, Y, labels)
    fig = Figure(size = (1200, 600))

    dim_X = size(X, 1)
    dim_Y = size(Y, 1)

    if dim_X >= 3
        ax_X = Axis3(fig[1, 1], title = "Przestrzeń X ($(dim_X)D)", azimuth = 0.22 * π)
        scatter!(ax_X, X[1, :], X[2, :], X[3, :], color = labels, colormap = :jet, markersize = 8)
    elseif dim_X == 2
        ax_X = Axis(fig[1, 1], title = "Przestrzeń X (2D)")
        scatter!(ax_X, X[1, :], X[2, :], color = labels, colormap = :jet, markersize = 8)
    else
        ax_X = Axis(fig[1, 1], title = "Przestrzeń X (1D)")
        # Dla 1D dodajemy wektor zer, by narysować punkty na płaskiej osi
        scatter!(ax_X, X[1, :], zeros(size(X, 2)), color = labels, colormap = :jet, markersize = 8)
    end

    if dim_Y >= 3
        ax_Y = Axis3(fig[1, 2], title = "Przestrzeń Y ($(dim_Y)D)")
        scatter!(ax_Y, Y[1, :], Y[2, :], Y[3, :], color = labels, colormap = :jet, markersize = 8)
    elseif dim_Y == 2
        ax_Y = Axis(fig[1, 2], title = "Przestrzeń Y (2D)")
        scatter!(ax_Y, Y[1, :], Y[2, :], color = labels, colormap = :jet, markersize = 8)
    else
        ax_Y = Axis(fig[1, 2], title = "Przestrzeń Y (1D)")
        scatter!(ax_Y, Y[1, :], zeros(size(Y, 2)), color = labels, colormap = :jet, markersize = 8)
    end

    return fig
end
