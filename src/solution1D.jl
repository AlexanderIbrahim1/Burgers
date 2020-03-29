# An example of the Burgers code for the propagation of Burgers equation in time
# plots u(x, t = t0) for a given t0
# studies the pointwise convergence order of a solution
# to Burgers equation as a function of resolution

using Burgers
using Printf
using PyPlot
include("Utils.jl")
using .Utils

function first_deriv(fun::FunArr{T, U})::FunArr{T, U} where {T, U}
	N     = fun.Ncells
	dudat = similar(fun.udat)
	
	# take first derivative of end-points
	dudat[1] = (fun.udat[2] - fun.udat[1]  ) / ufun.dx
	dudat[N] = (fun.udat[N] - fun.udat[N-1]) / ufun.dx
	
	# take the first derivative of all points in between
	for i in 2:(N-1)
		dudat[i] = (fun.udat[i+1] - fun.udat[i-1]) / (2*ufun.dx)
	end
	
	println(fun.Ncells, " ", fun.Nghost)
	FunArr{T, U}(dudat, fun.Ncells, fun.Nghost)
end

function second_deriv(fun::FunArr{T, U})::FunArr{T, U} where {T, U}
	N      = fun.Ncells
	d2udat = similar(fun.udat)
	
	# take second derivative of end-points
	d2udat[1] = (fun.udat[3] - 2.0*fun.udat[2]   + fun.udat[1])   / fun.dx^2
	d2udat[N] = (fun.udat[N] - 2.0*fun.udat[N-1] + fun.udat[N-2]) / fun.dx^2
	
	# take the second derivative of all points in between
	for i in 2:(N-1)
		d2udat[i] = (fun.udat[i+1] - 2*fun.udat[i] + fun.udat[i-1]) / fun.dx^2
	end
	
	println(fun.Ncells, " ", fun.Nghost)
	FunArr{T, U}(d2udat, fun.Ncells, fun.Nghost)
end

function plotsoln(fun::FunArr{T, U}, imagename::AbstractString, ylab::AbstractString) where {T, U}
	PyPlot.clf()
	xlabel("x")
	xlim(0, 1)
	xticks(range(0.0, 1.0, length = 10 + 1))
	
	ylabel(ylab)
	yticks(range(0.0, 0.7, length = 7 + 1))
	ylim(0, 0.7)
	plot(fun.xdat, fun.udat)
	savefig("images/" * imagename * ".png", dpi = 300)
end

# collect the command line arguments
t0, sc, imagename = commandlineparser(false)

# create the resolution
Ttot     = 2.0   # total duration
Nt_b     = 1024  # number of time divisions for lowest resolution
Ncells_b = 256   # number of position divisions for lowest resolution
Nghost   = 2     # padding on position divisions for boundary conditions
res = Resolution{Float64, Int64}(Ttot, Nt_b, Ncells_b, Nghost)
check_time(res, t0)

# find the evolved solution
ufun   = get_evolved_slice(get_sinewave, res, t0, sc)
dufun  = first_deriv(ufun)
d2ufun = second_deriv(ufun)

# plot u(x) and its derivatives and save the results
plotsoln(ufun  , imagename                 , "u(x)"  )
#plotsoln(dufun , imagename * "_firstderiv" , "du(x)" )
#plotsoln(d2ufun, imagename * "_secondderiv", "d2u(x)")
