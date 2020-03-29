# An example of the Burgers code for the
# propagation of Burgers equation in time
# studies the pointwise convergence order of a solution
# to Burgers equation as a function of resolution

using Burgers
using Printf
using PyPlot
include("Utils.jl")
using .Utils

function convergence_ratio(it::I, sc::SlopeCalc) where {T, I<:Int}
	# calculates the pointwise convergence ratio function
	# it     : the index for the time at which to slice u(x, t)
	#          must satisfy 0 <= t0 <= Ttot
	# sc     : the type of slope calculation (linear, minmod, constant)
	
	Ttot     = 2.0   # total duration
	Nt_b     = 512   # number of time divisions for lowest resolution
	Ncells_b = 128   # number of position divisions for lowest resolution
	Nghost   = 2     # padding on position divisions for boundary conditions
	
	if !(1 <= it <= Nt_b)
		error("The input time index must be between 1 and $Nt_base")
	end
	
	# The scaling factors
	scale1 = 8
	scale2 = 2*scale1
	scale3 = 2*scale2
	
	# Creating the Resolution instances
	res1 = Resolution{Float64, Int64}(Ttot, scale1*Nt_b, scale1*Ncells_b, Nghost)
	res2 = Resolution{Float64, Int64}(Ttot, scale2*Nt_b, scale2*Ncells_b, Nghost)
	res3 = Resolution{Float64, Int64}(Ttot, scale3*Nt_b, scale3*Ncells_b, Nghost)
	
	# Evaluating the solution to Burgers equation at t = t0
	ufun1 = get_evolved_slice(get_sinewave, res1, scale1*it, sc)
	ufun2 = get_evolved_slice(get_sinewave, res2, scale2*it, sc)
	ufun3 = get_evolved_slice(get_sinewave, res3, scale3*it, sc)
	
	# Calculating the absolute difference between solutions
	# of successive resolutions
	udiff12 = absolute_difference(ufun1, ufun2, ufun1.xdat)
	udiff23 = absolute_difference(ufun2, ufun3, ufun1.xdat)
	
	# stabilize the solutions, prevents the denominator from being too small
	eps = 1e-6
	uratio = (eps + udiff12)/(eps + udiff23)
	
	return uratio
end

function plot_convergence(cfun::FunArr{T, U}, imagename::AbstractString) where {T, U}
	# plots and saves the binary logarithm of the pointwise convergence function
	PyPlot.clf()
	
	xlabel("x")
	xlim(0, 1)
	xticks(range(0, 1, length = 10 + 1))

	ylabel("C(x)")
	ylim(0, 6)
	yticks(range(0, 6, length = 6 + 1))

	plot(cfun.xdat, log.(2, cfun.udat))
	savefig("images/" * imagename * ".png", dpi = 300)
end

# collect the command line arguments
# calculate and plot the convergence ratio
it, sc, imagename = commandlineparser(true)
convergence       = convergence_ratio(it, sc)
plot_convergence(convergence, imagename)
