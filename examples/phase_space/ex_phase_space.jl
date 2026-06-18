using DifferentialEquations
using GLMakie
using AtractorsQGP: set_publication_theme as theme_QGP
using LaTeXStrings

function ukos_oscylator(u, p, t)
    x, v = u
    c = p
    dx = v
    dv = -x + c * t
    return [dx, dv]
end

function generuj_wykresy_przestrzeni()

    tspan = (0.0, 10)
    warunki_poczatkowe = [
        [2.0, 0.0],
        [-2.0, 0.0],
        [1.0, 0.0],
        # [1.0, 1.0],
    ]
    theme_QGP()
    paleta = Makie.theme(:Palette).color[]


    fig1 = Figure(size = (850, 600))
    ax1_2d = Axis(fig1[1, 1], xlabel = L"x", ylabel = L"v", title = L"\text{Naiwna przestrzeń fazowa} \; 2D \; (c = 0)")

    fig2 = Figure(size = (850, 600))
    ax1_3d = Axis3(fig2[1, 1], xlabel = L"x", ylabel = L"v", zlabel = L"t", title = L"\text{Rozszerzona przestrzeń fazowa} \; 3D \; (c = 0)")

    for (i, u0) in enumerate(warunki_poczatkowe)
        prob = ODEProblem(ukos_oscylator, u0, tspan, 0.0)
        sol = solve(prob, Tsit5(), saveat = 0.01)
        x_pts = [u[1] for u in sol.u]
        v_pts = [u[2] for u in sol.u]

        lines!(ax1_2d, x_pts, v_pts, color = :gray, linewidth = 9)
        lines!(ax1_2d, x_pts, v_pts, color = paleta[5 * i], linewidth = 8)

        lines!(ax1_3d, x_pts, v_pts, sol.t, color = paleta[i * 5], linewidth = 7)
    end

    fig3 = Figure(size = (850, 600))
    ax2_2d = Axis(fig3[1, 1], xlabel = L"x", ylabel = L"v", title = L"\text{Naiwna przestrzeń fazowa} \; 2D \; (c = 0.5)")

    fig4 = Figure(size = (850, 600))
    ax2_3d = Axis3(fig4[1, 1], xlabel = L"x", ylabel = L"v", zlabel = L"t", title = L"\text{Rozszerzona przestrzeń fazowa} \; 3D \; (c = 0.5)")

    for (i, u0) in enumerate(warunki_poczatkowe)
        prob = ODEProblem(ukos_oscylator, u0, tspan, 0.5)
        sol = solve(prob, Tsit5(), saveat = 0.01)
        x_pts = [u[1] for u in sol.u]
        v_pts = [u[2] for u in sol.u]

        lines!(ax2_2d, x_pts, v_pts, color = :grey, linewidth = 10)
        lines!(ax2_2d, x_pts, v_pts, color = paleta[5 * i], linewidth = 8)

        lines!(ax2_3d, x_pts, v_pts, sol.t, color = paleta[i * 5], colormap = paleta, linewidth = 7)
    end

    save("plots/ex_phase_space_2d_c0.pdf", fig1)
    save("plots/ex_phase_space_3d_c0.pdf", fig2)
    save("plots/ex_phase_space_2d_c05.pdf", fig3)
    save("plots/ex_phase_space_3d_c05.pdf", fig4)

    return fig1, fig2, fig3, fig4
end
