# ==============================================================================
# Script: makie_extras.jl
# Contains additional/alternative visualization scripts using CairoMakie
# ==============================================================================

using CairoMakie
import LinearAlgebra: norm
using Random

"""
    plot_physical_trajectories()

Generates the "Physical Trajectories" plot containing a central glowing sphere (core)
and a distorted helical trajectory, using the corrected `cgrad` colormap.
Saves the output to `makie_physical_trajectories.png`.
"""
function plot_physical_trajectories()
    println("Generating: makie_physical_trajectories.png ...")
    set_theme!(theme_black())
    fig = Figure(size = (1000, 1000))
    ax = Axis3(fig[1, 1], aspect = :data, azimuth = 0.45, elevation = 0.5)
    hidespines!(ax)
    hidedecorations!(ax)

    central_pos = Point3f(0.0, 0.0, 0.0)

    # 1. Glow effect (multiple translucent spheres stacked together)
    # Using Sphere centered at origin so that meshscatter! correctly places it at central_pos
    for r in [0.5, 0.7, 1.0]
        meshscatter!(ax, [central_pos], marker=Sphere(Point3f(0), r),
                    color=(:white, 0.1), shading=true)
    end

    # Core (unsheared/unshaded core)
    meshscatter!(ax, [central_pos], marker=Sphere(Point3f(0), 0.2),
                color=:white, shading=false)

    # Helper function to generate distorted helix
    function generate_distorted_helix(t_range, radius, dist_coeff)
        X = [radius * cos(t) * exp(-dist_coeff*t) for t in t_range]
        Y = [radius * sin(t) * exp(-dist_coeff*t) for t in t_range]
        Z = [1.2 * t * exp(-0.01*t) for t in t_range] # Extra oscillation on Z
        
        return [Point3f(x, y, z) for (x,y,z) in zip(X, Y, Z)]
    end

    t_vals = range(-50, 50, length=2000)
    traj1 = generate_distorted_helix(t_vals, 3.0, 0.0)

    # Define custom colormap with transparency using RGBAf
    col_map = cgrad([:black, RGBAf(1.0, 1.0, 1.0, 0.3), :white])

    # Plot the helical trajectory
    lines!(ax, traj1, color = t_vals, colormap = col_map, linewidth=1.5, overdraw=true)

    save("makie_physical_trajectories.png", fig, px_per_unit = 2)
    println("Saved: makie_physical_trajectories.png")
    return fig
end

"""
    plot_waves_interferometry()

Generates a three-panel plot consisting of:
A. Point sources (using a custom `circles!` implementation)
B. 2D Interference field (using contour and wireframe)
C. Superimposed Signal Waves (Joyplot/waterfall)
Saves the output to `makie_waves_interferometry.svg`.
"""
function plot_waves_interferometry()
    println("Generating: makie_waves_interferometry.svg ...")
    set_theme!(theme_minimal())

    fig = Figure(size = (1200, 1000))
    ga = fig[1, 1] = GridLayout()

    # --- 1. Point Sources (Top-left) ---
    ax_dots = Axis(ga[1, 1], title="A. Punktowe źródła")
    hidespines!(ax_dots)
    hidedecorations!(ax_dots)

    # Custom circles! function implementation since it isn't in standard Makie
    function local_circles!(ax, xs, ys; radii, color, strobe=true)
        for (x, y) in zip(xs, ys)
            for r in radii
                arc!(ax, Point2f(x, y), r, 0.0, 2π, color = color, linewidth = 1.0)
            end
        end
    end

    n_dots = 4
    Random.seed!(42) # Set seed for reproducibility of random radii
    for i in 1:n_dots, j in 1:n_dots
        r = 0.1 * abs(randn())
        local_circles!(ax_dots, [j], [i], radii=[r, r*1.5, r*2.0], 
                      color=(:black, 0.3), strobe=true)
        scatter!(ax_dots, [j], [i], color=:black, markersize=10)
    end

    # --- 2. 2D Interference Field (Top-right) ---
    ax_wave = Axis3(ga[1, 2], title="B. Pole Interferencyjne")
    hidespines!(ax_wave)
    hidedecorations!(ax_wave)

    xs = range(-5, 5, length=100)
    ys = range(-5, 5, length=100)
    zs = [sin(norm([x, y])) / (0.1 + norm([x, y])^2) for x in xs, y in ys]
    
    contour!(ax_wave, xs, ys, zs, color=(:black, 0.4), linewidth=0.5, levels=30)
    wireframe!(ax_wave, xs, ys, zs, color=(:black, 0.1), linewidth=0.2)

    # --- 3. Superimposed Signal Waves / Joyplot (Bottom) ---
    ax_sig = Axis(ga[2, 1:2], title="C. Nakładanie się Sygnałów (Joyplot)")
    hidespines!(ax_sig)
    hidedecorations!(ax_sig)

    t = range(0, 10, length=1000)
    num_sigs = 15

    for i in 1:num_sigs
        envelope = exp.(-(t .- 5 .- (i * 0.2)).^2 ./ 4)
        sig = 0.5 .* envelope .* sin.((1 + 0.1 * i) .* t .+ i * 0.5)
        
        # Add a bit of noise (optional, kept for organic feeling as per original code)
        noise = 0.05 .* cumsum(randn(length(t)))
        smoothed_noise = exp.(-t.^2/2) .* noise
        
        y_offset = num_sigs - i
        
        # Draw the lines
        lines!(ax_sig, t, sig .+ y_offset, color=(:black, 0.6), linewidth=1.1)
    end

    save("makie_waves_interferometry.svg", fig)
    println("Saved: makie_waves_interferometry.svg")
    return fig
end

"""
    plot_torus_minimal()

Generates a minimal, elegant wireframe of a torus on a dark background.
Saves the output to `makie_torus_minimal.png`.
"""
function plot_torus_minimal()
    println("Generating: makie_torus_minimal.png ...")
    set_theme!(theme_black())

    # Function generating points for Torus
    function generate_torus(R, r; resolution=60)
        U = range(0, 2π, length=resolution)
        V = range(0, 2π, length=resolution)
        X = [ (R + r*cos(v))*cos(u) for u in U, v in V ]
        Y = [ (R + r*cos(v))*sin(u) for u in U, v in V ]
        Z = [ r*sin(v) for u in U, v in V ]
        return X, Y, Z
    end

    fig = Figure(size = (1000, 1000))
    ax = Axis3(fig[1, 1], aspect = :data, azimuth = 0.45, elevation = 0.5)
    hidespines!(ax)
    hidedecorations!(ax)

    R_val = 3.0
    r_val = 1.0
    X, Y, Z = generate_torus(R_val, r_val; resolution=70)

    # Main wireframe
    wireframe!(ax, X, Y, Z, color = (:white, 0.3), linewidth = 0.7, overdraw=true)

    # Additional lighter outer grid for depth effect
    X_out, Y_out, Z_out = generate_torus(R_val, r_val; resolution=10)
    wireframe!(ax, X_out, Y_out, Z_out, color = (:white, 0.1), linewidth = 1.2)

    save("makie_torus_minimal.png", fig, px_per_unit = 2)
    println("Saved: makie_torus_minimal.png")
    return fig
end

# If run directly as a script:
if abspath(PROGRAM_FILE) == @__FILE__
    plot_physical_trajectories()
    plot_waves_interferometry()
    plot_torus_minimal()
end
