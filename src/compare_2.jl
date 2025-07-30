using HydroAttractors
using Plots

u0 = [2.4228426395939087, 5.5]
u1 = [2.4228426395939087, 2.5]

sett = SimSettings(tspan=(0.1, 1))

sol0 = evol(u0, sett)
sol1 = evol(u1, sett)

using Plots

p = plot(sol0, idxs=2, xlabel="t", ylabel="A")
plot!(sol1, idxs=2, xlabel="t", ylabel="A")


display(p)
println("Press Enter ...")
readline()


 
