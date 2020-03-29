# An example of the Burgers code for the
# propagation of Burgers equation in time
# Creates a colormap showing u(x, t) of Burgers equation
# initialized with the sine function used in lecture

using Burgers
using Printf
using PyPlot
include("Utils.jl")
using .Utils

Ttot   = 2.0   # total duration
Nt     = 1024  # number of time divisions
Ncells = 256   # number of position divisions
Nghost = 2     # padding on position divisions for boundary conditions

res   = Resolution{Float64, Int64}(Ttot, Nt, Ncells, Nghost)
sc    = minmod::SlopeCalc
u2fun = get_evolved(get_sinewave, res, sc)

# u(x, t) can be plotted as a colour map
xlabel("t")
ylabel("x")
pcolormesh(u2fun.tdat, u2fun.xdat, u2fun.u2dat)
colorbar()
savefig("images/colormap2.png", dpi = 300)
