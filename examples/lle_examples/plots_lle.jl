
function plot_examples_lle(X, Y, labels)
    fig = Figure(size = (1200, 600))
    sizeX = size(X)
    ax_3d = Axis3(fig[1, 1], title = "Oryginalna rozmaitość X $sizeX", azimuth = 0.22 * π)
    scatter!(ax_3d, X[1, :], X[2, :], X[3, :], color = labels, colormap = :jet, markersize = 8)
    # println(Y)
    ax_2d = Axis3(fig[1, 2], title = "Zredukowana przestrzeń Y $(size(Y))", azimuth = 0.22 * π)
    scatter!(ax_2d, Y[1, :], Y[2, :], Y[3,:],color = labels, colormap = :jet, markersize = 8)
    
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






function set_publication_theme_large()
    set_theme!(
        Theme(
            font = "Libertinus", 
            fontsize = 35,           
            figure_padding = 30,
            
            Axis = (
                titlesize = 50,
                xlabelsize = 35,
                ylabelsize = 35,
                xticklabelsize = 28,
                yticklabelsize = 28,
                backgroundcolor = RGBf(1.0, 1.0, 1.0),
                xgridstyle = :dash,
                ygridstyle = :dash,
                xgridcolor = RGBAf(0.85, 0.85, 0.85, 0.65),
                ygridcolor = RGBAf(0.85, 0.85, 0.85, 0.65),
                spinewidth = 2.0,     
                xtickwidth = 2.0,
                ytickwidth = 2.0,
                xgridwidth = 1.5,
                ygridwidth = 1.5,
                
                xticksize = 10,       
                yticksize = 10,
                xtickalign = 1.0,
                ytickalign = 1.0,
                topspinevisible = true,
                rightspinevisible = true,
            ),
            
            Legend = (
                titlesize = 28,
                labelsize = 24,
                framevisible = true,
                framewidth = 1.5,     # Grubsza ramka legendy
                framecolor = RGBAf(0.8, 0.8, 0.8, 1.0),
                backgroundcolor = RGBAf(1.0, 1.0, 1.0, 0.85),
                position = :rt,
                padding = (12.0, 12.0, 12.0, 12.0), # Więcej oddechu wewnątrz
            ),
            Palette = (color = Makie.wong_colors(),),
            Lines = (
                linewidth = 3.0,
            ),
            Scatter = (
                markersize = 12,
            )
        )
    )
end

