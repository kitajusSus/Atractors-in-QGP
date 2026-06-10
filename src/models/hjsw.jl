"""
Konfiguracja modelu HJSW (Heller-Janik-Spaliński-Witaszczyk).
Zawiera klasyczne parametry transportu oraz stałe tłumienia modów QNM.
"""
struct HJSWModel{T <: Real} <: AbstractHydroModel
    params::HydroParams{T}
    omega_R::T
    omega_I::T
end

"""
    HJSWModel(; eta_over_s, tau_pi, lambda1, omega_R, omega_I)

Tworzy model HJSW z domyślnymi wartościami z artykułu dla plazmy N=4 SYM.
"""
function HJSWModel(;
        eta_over_s::Real = 1 / (4 * π),
        tau_pi::Real = 1 / (2 * π),
        lambda1::Real = 0.0,
        omega_R::Real = 9.8,
        omega_I::Real = 8.629,
    )
    T = promote_type(typeof(eta_over_s), typeof(tau_pi), typeof(lambda1), typeof(omega_R), typeof(omega_I))
    params = HydroParams(T(eta_over_s), T(tau_pi), T(lambda1))
    return HJSWModel(params, T(omega_R), T(omega_I))
end
