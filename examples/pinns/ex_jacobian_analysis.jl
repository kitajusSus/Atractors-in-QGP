"""
    ex_jacobian_analysis.jl

Example: Jacobian-based dimensionality reduction via parametric PINN.

Demonstrates the full pipeline:
1. Train a parametric PINN for Bjorken flow (BRSSS/MIS)
2. Compute the Jacobian  J(τ) = ∂(T̂, Â)/∂(T₀, A₀)  via ForwardDiff
3. Track d_eff(τ) = (Σσᵢ²)²/Σσᵢ⁴ over proper time
4. Detect hydrodynamisation time τ_hyd where d_eff → 1

For multidimensional data:
- `pinn_jacobian_full` gives the 2×4 Jacobian including η/s and τ_π
- `pinn_jacobian_scan` returns raw J matrices for all (τ, IC) pairs
"""

using AttractorsQGP
using GLMakie
using Printf
using Statistics

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Configuration and training
# ─────────────────────────────────────────────────────────────────────────────

println("="^65)
println("  PINN Jacobian Analysis — QGP Dimensionality Reduction")
println("="^65)

model = BRSSSModel()   # BRSSS hydrodynamics
config = PINNConfig(
    hidden_layers = 5,
    hidden_size = 128,
    τ_range = (0.15, 1.5),
    T_range = (0.5, 15.0),
    A_range = (-15.0, 25.0),
    η_range = (0.05, 1.0),
    τπ_range = (0.05, 1.0),
)

println("\n[1/4] Training PINN …")
result = train_pinn(
    model, config;
    n_epochs = 6000,
    batch_size_colloc = 512,
    batch_size_ic = 128,
    learning_rate = 1.0e-3,
    λ_physics = 1.0,
    λ_ic = 10.0,
    verbose = true,
    log_every = 200,
)

# ─────────────────────────────────────────────────────────────────────────────
# 2.  Build IC ensembles
# ─────────────────────────────────────────────────────────────────────────────

println("\n[2/4] Building IC ensembles …")

# Full ensemble — samples η/s and τ_π too
ics_full = sample_ic_ensemble(result, 200; seed = 7)

# Fixed transport — only state ICs vary (useful for cleaner comparison)
ics_fixed = fixed_transport_ic_ensemble(
    result, 200;
    η = model.params.eta_over_s,
    τπ = model.params.tau_pi,
    seed = 42,
)

taus = collect(0.2:0.02:1.4)

# ─────────────────────────────────────────────────────────────────────────────
# 3.  d_eff scan  (2×2 Jacobian w.r.t. state ICs only)
# ─────────────────────────────────────────────────────────────────────────────

println("\n[3/4] Running Jacobian d_eff scan (state ICs) …")
scan_state = pinn_deff_scan(result, taus, ics_fixed; full = false, verbose = true)

τ_hyd_state = pinn_hydrodynamisation_time(scan_state; threshold = 1.05)
@printf("  → τ_hyd (state)  = %.3f fm/c\n", τ_hyd_state)

# ─────────────────────────────────────────────────────────────────────────────
# 4.  d_eff scan  (2×4 Jacobian w.r.t. all 4 input parameters)
# ─────────────────────────────────────────────────────────────────────────────

println("\n[4/4] Running Jacobian d_eff scan (full: state + transport) …")
scan_full = pinn_deff_scan(result, taus, ics_full; full = true, verbose = true)

τ_hyd_full = pinn_hydrodynamisation_time(scan_full; threshold = 1.05)
@printf("  → τ_hyd (full)   = %.3f fm/c\n", τ_hyd_full)

# ─────────────────────────────────────────────────────────────────────────────
# Optional: raw Jacobian scan for geometry
# ─────────────────────────────────────────────────────────────────────────────

# Only compute for a small sub-ensemble (it returns all J matrices)
ics_small = ics_fixed[1:20]
jac_scan = pinn_jacobian_scan(
    result, taus[1:10:end], ics_small;
    full = false, verbose = false
)

# Largest singular value at τ=taus[1] for first IC
σ_example = jac_scan.svd_values[1, 1]
@printf("  σ₁=%.4f  σ₂=%.4f  (τ=%.2f, IC #1)\n", σ_example[1], σ_example[2], taus[1])

# ─────────────────────────────────────────────────────────────────────────────
# Visualisation
# ─────────────────────────────────────────────────────────────────────────────

println("\nGenerating figures …")

fig = Figure(size = (1200, 800), fontsize = 14)

# ── Panel 1: d_eff(τ) — state Jacobian ──────────────────────────────────────
ax1 = Axis(
    fig[1, 1],
    title = "Effective Dimension d_eff(τ) — State Jacobian ∂(T,A)/∂(T₀,A₀)",
    xlabel = "Proper time τ [fm/c]",
    ylabel = "d_eff(τ)",
)
band!(
    ax1, scan_state.taus,
    scan_state.d_eff_mean .- scan_state.d_eff_std,
    scan_state.d_eff_mean .+ scan_state.d_eff_std,
    color = (:royalblue, 0.25)
)
lines!(
    ax1, scan_state.taus, scan_state.d_eff_mean,
    color = :royalblue, linewidth = 3, label = "d_eff (state)"
)
hlines!(ax1, [2.0], color = :gray, linestyle = :dash, label = "2D (far from equilibrium)")
hlines!(ax1, [1.0], color = :crimson, linestyle = :dash, label = "1D attractor")
if isfinite(τ_hyd_state)
    vlines!(
        ax1, [τ_hyd_state], color = :green, linestyle = :dashdot, linewidth = 2,
        label = @sprintf("τ_hyd = %.2f fm/c", τ_hyd_state)
    )
end
axislegend(ax1; position = :rt)
ylims!(ax1, 0.8, 2.3)

# ── Panel 2: d_eff(τ) — full Jacobian ───────────────────────────────────────
ax2 = Axis(
    fig[1, 2],
    title = "Effective Dimension d_eff(τ) — Full Jacobian ∂(T,A)/∂(T₀,A₀,η,τπ)",
    xlabel = "Proper time τ [fm/c]",
    ylabel = "d_eff(τ)",
)
band!(
    ax2, scan_full.taus,
    scan_full.d_eff_mean .- scan_full.d_eff_std,
    scan_full.d_eff_mean .+ scan_full.d_eff_std,
    color = (:darkorange, 0.25)
)
lines!(
    ax2, scan_full.taus, scan_full.d_eff_mean,
    color = :darkorange, linewidth = 3, label = "d_eff (full)"
)
hlines!(ax2, [2.0], color = :gray, linestyle = :dash, label = "2D limit")
hlines!(ax2, [1.0], color = :crimson, linestyle = :dash, label = "1D attractor")
if isfinite(τ_hyd_full)
    vlines!(
        ax2, [τ_hyd_full], color = :green, linestyle = :dashdot, linewidth = 2,
        label = @sprintf("τ_hyd = %.2f fm/c", τ_hyd_full)
    )
end
axislegend(ax2; position = :rt)
ylims!(ax2, 0.8, 2.3)

# ── Panel 3: Singular values σ₁(τ), σ₂(τ) ──────────────────────────────────
ax3 = Axis(
    fig[2, 1],
    title = "Singular Values of J(τ) — State Jacobian",
    xlabel = "Proper time τ [fm/c]",
    ylabel = "σᵢ(τ)  [log scale]",
    yscale = log10,
)
lines!(
    ax3, scan_state.taus, scan_state.sigma1_mean,
    color = :royalblue, linewidth = 3, label = "σ₁(τ)"
)
lines!(
    ax3, scan_state.taus, max.(scan_state.sigma2_mean, fill(1.0e-12, length(taus))),
    color = :crimson, linewidth = 3, linestyle = :dash, label = "σ₂(τ)"
)
axislegend(ax3; position = :rt)

# ── Panel 4: σ₂/σ₁ ratio ────────────────────────────────────────────────────
ax4 = Axis(
    fig[2, 2],
    title = "Singular Value Ratio σ₂/σ₁ — Attractor Indicator",
    xlabel = "Proper time τ [fm/c]",
    ylabel = "σ₂(τ) / σ₁(τ)",
)
lines!(
    ax4, scan_state.taus, scan_state.sigma_ratio,
    color = :purple, linewidth = 3, label = "σ₂/σ₁ (state)"
)
lines!(
    ax4, scan_full.taus, scan_full.sigma_ratio,
    color = :darkorange, linewidth = 3, linestyle = :dash, label = "σ₂/σ₁ (full)"
)
hlines!(ax4, [0.0], color = :black, linestyle = :dot)
axislegend(ax4; position = :rt)

# Save
outpath = joinpath(@__DIR__, "jacobian_deff_analysis.png")
save(outpath, fig)
println("  Saved: ", outpath)

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "="^65)
println("  SUMMARY")
println("="^65)
@printf("  Model            : %s\n", typeof(model))
@printf("  η/s              : %.3f\n", model.params.eta_over_s)
@printf("  τ_π              : %.3f fm/c\n", model.params.tau_pi)
@printf("  IC ensemble      : %d points\n", length(ics_fixed))
@printf("  τ_hyd (state J)  : %.3f fm/c\n", τ_hyd_state)
@printf("  τ_hyd (full J)   : %.3f fm/c\n", τ_hyd_full)
@printf("  d_eff at τ=%.2f : %.4f\n", taus[1], scan_state.d_eff_mean[1])
@printf("  d_eff at τ=%.2f : %.4f\n", taus[end], scan_state.d_eff_mean[end])
println("="^65)
