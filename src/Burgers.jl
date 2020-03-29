module Burgers

using DataInterpolations

export FunArr
export FunArr2D
export evolve
export get_timeslice
export strip_boundary_conditions
export integrate_difference1D
export absolute_difference
export array_smoother
export SlopeCalc
export minmod
export linear
export constant

@enum SlopeCalc minmod linear constant

mutable struct FunArr2D{T, U}
	xmin::T
	xmax::T
	dx::T
	Nx::Int64
	
	tmin::T
	tmax::T
	dt::T
	Nt::Int64
	
	xdat::Vector{T}
	tdat::Vector{T}
	u2dat::Array{U, 2}
	
	function FunArr2D{T, U}(_xdat, t_duration, _Nt) where {T, U}
		# test for appropriate sizes
		if length(_xdat) < 2
			error("The first dimension of the array must be at least 2.")
		end
		if _Nt < 2
			error("The second dimension of the array must be at least 2.")
		end
		
		# the duration of elapsed time must be positive
		if t_duration <= 0.0
			error("The duration of time that passes must be positive.")
		end
		
		_xmin = _xdat[1]
		_xmax = _xdat[end]
		_dx   = _xdat[2] - _xdat[1]
		_Nx   = length(_xdat)
		_tmin = 0.0
		_tmax = t_duration
		_dt = (_tmax - _tmin)/(_Nt - 1)
		_tdat = range(0.0, stop = t_duration, length = _Nt)
		_u2dat = Array{U, 2}(undef, _Nx, _Nt)
		
		new(_xmin, _xmax, _dx, _Nx, _tmin, _tmax, _dt, _Nt, _xdat, _tdat, _u2dat)
	end
end

mutable struct FunArr{T, U}
	xmin::T            # x-domain minimum
	xmax::T            # x-domain maximum
	xmin_g::T          # x-domain minimum with ghost
	xmax_g::T          # x-domain maximum with ghost
	Ncells::Int64      # size of the array held, excluding ghost cells
	Nghost::Int64      # number of ghost elements at EACH END of the array
	Nsize::Int64       # the total size of the array, with ghost cells included
	dx::T              # step size of z-data
	xdat::Vector{T}    # domain data
	udat::Vector{U}    # range data
	
	function FunArr{T, U}(func::Function, _Ncells, _Nghost) where {T, U}
		# func : function at which the xdat array is evaluated
		#      : must have the signature func(::Vector{T})::Vector{T}
		_xmin   = 0.0
		_xmax   = 1.0
		_dx     = (_xmax - _xmin)/(_Ncells)
		_xmin_g = _xmin + T(_dx*(0.5 - _Nghost))
		_xmax_g = _xmax + T(_dx*(_Nghost - 0.5))
		_Nsize  = _Ncells + 2*_Nghost
		_xdat   = range(_xmin_g, stop = _xmax_g, length = _Nsize)
		_udat   = func(_xdat)
		new(_xmin, _xmax, _xmin_g, _xmax_g, _Ncells, _Nghost, _Nsize, _dx, _xdat, _udat)
	end
	
	function FunArr{T, U}(_udat, _Ncells, _Nghost) where {T, U}
		_xmin   = 0.0
		_xmax   = 1.0
		_dx     = (_xmax - _xmin)/(_Ncells)
		_xmin_g = _xmin + _dx*(0.5 - T(_Nghost))
		_xmax_g = _xmax + _dx*(T(_Nghost) - 0.5)
		_Nsize  = _Ncells + 2*_Nghost
		_xdat   = range(_xmin_g, stop = _xmax_g, length = _Nsize)
		new(_xmin, _xmax, _xmin_g, _xmax_g, _Ncells, _Nghost, _Nsize, _dx, _xdat, _udat)
	end
end

function Base.:+(a::FunArr{T, U}, b::FunArr{T, U})::FunArr{T, U} where {T, U}
	# addition for two FunArr{T, U} instances
	if (a.Ncells != b.Ncells) || (a.Nghost != b.Nghost)
		error("Two FunArr instances must hold vectors of the same size for addition.")
	end
	
	_udat = a.udat + b.udat
	FunArr{T, U}(_udat, b.Ncells, b.Nghost)
end

function Base.:+(a::Number, b::FunArr{T, U})::FunArr{T, U} where {T, U}
	# addition for two FunArr{T, U} instances
	_udat = a .+ b.udat
	FunArr{T, U}(_udat, b.Ncells, b.Nghost)
end

function Base.:-(a::FunArr{T, U}, b::FunArr{T, U})::FunArr{T, U} where {T, U}
	# subtraction for two FunArr{T, U} instances
	if (a.Ncells != b.Ncells) || (a.Nghost != b.Nghost)
		error("Two FunArr instances must hold vectors of the same size for subtraction.")
	end
	
	_udat = a.udat - b.udat
	FunArr{T, U}(_udat, b.Ncells, b.Nghost)
end

function Base.:*(a::FunArr{T, U}, b::Number)::FunArr{T, U} where {T, U}
	# multiplication of a FunArr{T, U} instance by a Number
	_udat = a.udat * b
	FunArr{T, U}(_udat, a.Ncells, a.Nghost)
end

function Base.:*(a::Number, b::FunArr{T, U})::FunArr{T, U} where {T, U}
	# multiplication of a FunArr{T, U} instance by a Number
	_udat = a * b.udat
	FunArr{T, U}(_udat, b.Ncells, b.Nghost)
end

function Base.:*(a::FunArr{T, U}, b::FunArr{T, U})::FunArr{T, U} where {T, U}
	# multiplication of a FunArr{T, U} instance by a Number
	_udat = a.udat .* b.udat
	FunArr{T, U}(_udat, b.Ncells, b.Nghost)
end

function Base.:/(a::FunArr{T, U}, b::FunArr{T, U})::FunArr{T, U} where {T, U}
	# multiplication of a FunArr{T, U} instance by a Number
	_udat = a.udat ./ b.udat
	FunArr{T, U}(_udat, b.Ncells, b.Nghost)
end

function Base.:/(a::FunArr{T, U}, b::Number)::FunArr{T, U} where {T, U}
	# multiplication of a FunArr{T, U} instance by a Number
	_udat = a.udat / b
	FunArr{T, U}(_udat, a.Ncells, a.Nghost)
end

###############################################################################
###############################################################################
###############################################################################

function throw_hasNaN(uarr::Vector{T}, msg) where {T}
	if sum(isnan.(uarr)) > 0
		println(uarr)
		error(msg)
	end
end

function physFlux(ufun::Vector{T})::Vector{T} where {T}
	# returns the value of F in the PDE:
	# (dt)u + (dx)F = 0
	return 0.5 * ufun .* ufun
end


function flux_naive(uarrP::Vector{T}, uarrM::Vector{T})::Vector{T} where {T}
	# returns the uncorrected average flux function at a boundary between cells
	return 0.5*(physFlux(uarrP) + physFlux(uarrM))
end

function flux_roe(uarrP::Vector{T}, uarrM::Vector{T})::Vector{T} where {T}
	# returns the average flux function at a boundary between cells
	# with the Roe's method correction, using an arithmetic average
	Fluxfun_naive = 0.5*(physFlux(uarrP) + physFlux(uarrM))
	uavgfun       = 0.5*abs.(uarrP + uarrM)
	ustepfun      = uarrM - uarrP
	return Fluxfun_naive - 0.5 * uavgfun .* ustepfun
end


function get_minmod_plusminus(uarr::Vector{T}) where {T}
	# uses the min-mod function to calculate the interface
	N = length(uarr)
	
	throw_hasNaN(uarr, "uarr")
	
	# the forward difference slopes
	sF      = Vector{T}(undef, N)
	sF[end] = uarr[1] - uarr[end]
	for i in 1:(N-1)
		sF[i] = uarr[i+1] - uarr[i]
	end

	# the backward difference slopes
	sB      = Vector{T}(undef, N)
	sB[1]   = uarr[1] - uarr[end]
	for i in 2:N
		sB[i] = uarr[i] - uarr[i-1]
	end
	
	# minmod-applied slope
	s = Vector{T}(undef, N)
	for i in 1:N
		if sF[i]*sB[i] < 0.0
			s[i] = 0.0
		elseif abs(sF[i]) < abs(sB[i])
			s[i] = sF[i]
		elseif abs(sB[i]) <= abs(sF[i])
			s[i] = sB[i]
		else
			error("This condition shouldn't happen, but just in case.")
		end
	end
	
	uarr_plus  = uarr + 0.5*s
	uarr_minus = uarr - 0.5*s
	
	return (uarr_plus, uarr_minus)
end

function get_linear_plusminus(uarr::Vector{T}) where {T}
	# "linear interpolation to interface"
	uarr_plus  = copy(uarr)
	uarr_minus = copy(uarr)
	uarr_plus[1:end-1] = 0.5*(uarr[2:end] + uarr[1:end-1])
	uarr_minus[2:end]  = 0.5*(uarr[2:end] + uarr[1:end-1])
	return (uarr_plus, uarr_minus)
end

function get_constant_plusminus(uarr::Vector{T}) where {T}
	# "linear interpolation to interface"
	uarr_plus  = copy(uarr)
	uarr_minus = copy(uarr)
	return (uarr_plus, uarr_minus)
end

function get_plusminus(uarr::Vector{T}, sc::SlopeCalc) where {T}
	# choose which method of calculating uarr_plus and uarr_minus to use
	if     sc == minmod
		return get_minmod_plusminus(uarr)
	elseif sc == linear
		return get_linear_plusminus(uarr)
	elseif sc == constant
		return get_constant_plusminus(uarr)
	else
		error("Invalid SlopeCalc instance used.")
	end
end


function udot(ufun::FunArr{T, U}, dx::T, sc::SlopeCalc)::FunArr{T, U} where {T, U}
	# used in the calculation of each step of RK2
	# finding ufun_plus and ufun_minus
	uarr_plus, uarr_minus = get_plusminus(ufun.udat, sc)
	
	Farr = flux_roe(uarr_plus[1:end-1], uarr_minus[2:end])
	Farr = vcat([0], Farr, [0])
	udot_arr = (Farr[1:end-1] - Farr[2:end])/dx
	udot_fun = 1.0*ufun
	udot_fun.udat = udot_arr
	
	throw_hasNaN(udot_arr, "udot_arr")
	
	return udot_fun
end

function step_RK2(ufun::FunArr{T, U}, dt::T, dx::T, sc::SlopeCalc)::FunArr{T, U} where {T, U}
	# performs a Runge-Kutta 2 step in time
	
	# first half-step of RK2
	ufun_half = ufun + 0.5*dt*udot(ufun, dx, sc)
	ufun_half = apply_boundary_conditions(ufun_half)
	throw_hasNaN(ufun_half.udat, "ufun_half.udat")
	
	# second half-step of RK2
	ufun_full = ufun + dt*udot(ufun_half, dx, sc)
	ufun_full = apply_boundary_conditions(ufun_full)
	throw_hasNaN(ufun_full.udat, "ufun_full.udat")
	
	return ufun_full
end


function strip_boundary_conditions(ufun::FunArr{T, U})::Vector{T} where {T, U}
	# removes the ghost cells from the ends of a FunArr instance
	Nghost = ufun.Nghost
	uarr   = ufun.udat[1+Nghost:end-Nghost]
	return uarr
end

function apply_boundary_conditions(ufun::FunArr{T, U})::FunArr{T, U} where {T, U}
	# copies ghost cells from one end to the other
	# assuming periodic boundary conditions
	uarr_bc = copy(ufun.udat)
	Nghost  = ufun.Nghost
	uarr_bc[1:Nghost]         = ufun.udat[end-2*Nghost+1:end-Nghost]
	uarr_bc[end-Nghost+1:end] = ufun.udat[Nghost+1:2*Nghost]
	
	ufun_bc = 1.0*ufun  # copy instead of point? IDK
	ufun_bc.udat[1:Nghost]         = ufun.udat[end-2*Nghost+1:end-Nghost]
	ufun_bc.udat[end-Nghost+1:end] = ufun.udat[Nghost+1:2*Nghost]
	return ufun_bc
end

function evolve(ufun::FunArr{T, U}, Nt::I, duration::T, sc::SlopeCalc, itstop::I) where {T, U, I<:Int}
	# ufun  : initial positions at t = 0
	# Nt    : total number of time steps
	# duration : total amount of time elapsed
	# returns : FunArr2D instance for the u(x, t) grid, and an array of energies
	
	dx = ufun.dx          # position step
	dt = duration/(Nt-1)  # time step
	
	Nghost = ufun.Nghost
	xarr  =  ufun.xdat[1+Nghost:end-Nghost]
	u2fun = FunArr2D{T, U}(xarr, duration, Nt)
	
	u2fun.u2dat[:,1] = strip_boundary_conditions(ufun)
	
	# propagate forward in time (need to make)
	ufun_next = ufun
	for it in 2:Nt
		ufun_next = step_RK2(ufun_next, dt, dx, sc)
		u2fun.u2dat[:,it] = strip_boundary_conditions(ufun_next)
		if it == itstop
			break
		end
	end
	
	return u2fun
end

function nearest_lower_index(N::I, xmin::T, xmax::T, x::T) where {T, I<:Int}
	# takes care of the rand() = 0.0 corner case
	if x == xmin
		x += eps(1.0)
	elseif x == xmax
		x -= eps(1.0)
	end
	if !(xmax > x > xmin)
		error("x = $x")
	end
	ratio = (x - xmin)/(xmax - xmin)
	lower_index = ceil(typeof(N), (N-1)*ratio)
	return lower_index
end

function get_timeslice(u2fun::FunArr2D{T, U}, t::T, Ncells::I, Nghost::I)::FunArr{T, U} where {T, U, I<:Int}
	# accepts a u2fun containing a time-evolved grid, and a time-step
	# returns the function u(x) evaluated at the timeslice it
	
	if !(u2fun.tmin <= t <= u2fun.tmax)
		error("Time out of bounds for this FunArr2D instance.")
	end
	
	NtL    = nearest_lower_index(u2fun.Nt, u2fun.tmin, u2fun.tmax, t)
	uarr1  = copy(u2fun.u2dat[:,NtL])
	uarr2  = copy(u2fun.u2dat[:,NtL+1])
	tratio = (t - u2fun.tdat[NtL])/u2fun.dt
	udat   = uarr1 .+ (uarr2 .- uarr1)*tratio
	ufun   = FunArr{T, U}(udat, Ncells, 0)
	
	return ufun
end

###############################################################################
# FUNCTIONS FOR CALCULATING L2 DIFFERENCE BETWEEN FUNCTIONS
###############################################################################

#function absolute_difference(ufun1::FunArr{T, U}, ufun2::FunArr{T, U}, Ndiv::I)::FunArr{T, U} where {T, U, I<:Int}
function absolute_difference(ufun1::FunArr{T, U}, ufun2::FunArr{T, U}, xdata::Vector{T})::FunArr{T, U} where {T, U, I<:Int}
	# returns the absolute difference between the linearly-extrapolated function ufun1 and ufun2
	@assert ufun1.xmin == ufun2.xmin
	@assert ufun1.xmax == ufun2.xmax
	
	xmin = ufun1.xmin
	xmax = ufun1.xmax
	Nghost = ufun1.Nghost
	
	interp1 = CubicSpline(ufun1.udat, ufun1.xdat)
	interp2 = CubicSpline(ufun2.udat, ufun2.xdat)
	
	uarr_diff = Vector{T}(undef, length(xdata))
	for (i, x) in enumerate(xdata)
		u1 = interp1(x)
		u2 = interp2(x)
		if (u1 == nothing) || (u2 == nothing)
			uarr_diff[i] = 0.0
		else
			uarr_diff[i] = abs(u2 - u1)
		end
	end
	
	return FunArr{T, U}(uarr_diff, length(xdata), 0)
end

function array_smoother(ufun::FunArr{T, U}, Nsmooth::I)::FunArr{T, U} where {T, U, I<:Int}
	# returns an array ufun_smooth where ufun_smooth[i] = avg(ufun, (i-Nsmooth):(i+Nsmooth))
	
	Nsize = ufun.Nsize
	uarr_smooth = Vector{T}(undef, Nsize)
	
	for i in 1:Nsize
		if ((i - Nsmooth) < 1) || ((i + Nsmooth) > Nsize)
			uarr_smooth[i] = ufun.udat[i]
			continue
		end
		
		avg = sum(ufun.udat[i-Nsmooth:i+Nsmooth])/(1 + 2*Nsmooth)
		uarr_smooth[i] = avg
	end
	
	ufun_smooth = FunArr{T, U}(uarr_smooth, Nsize, 0)
end

end # module

#function interp_regular1D(ufun::FunArr{T, U}, x::T) where {T, U}
#	# performs linear interpolation on a function discretized on a 1D regular grid
#	NxL = nearest_lower_index(ufun.Nsize, ufun.xmin, ufun.xmax, x)
#	x1 = ufun.xmin + (NxL-1)*ufun.dx
#	x2 = ufun.xmin + NxL*ufun.dx
#	u1 = ufun.udat[NxL]
#	u2 = ufun.udat[NxL+1]
#	slope = (u2 - u1)/ufun.dx
#	
#	interp_value = u1 + slope*(x-x1)
#	return interp_value
#end
#
#function integrate_difference1D(ufun1::FunArr{T, U}, ufun2::FunArr{T, U}, Ndiv::I) where {T, U, I<:Int}
#	@assert ufun1.xmin == ufun2.xmin
#	@assert ufun1.xmax == ufun2.xmax
#	
#	xmin = ufun1.xmin
#	xmax = ufun1.xmax
#	dx   = (xmax-xmin)/T(Ndiv-1)
#	
#	total = 0.0
#	for i in 1:(Ndiv-1)
#		x = (i-0.5)*dx
#		term_fun1 = interp_regular1D(ufun1, x)
#		term_fun2 = interp_regular1D(ufun2, x)
#		total    += dx * (term_fun1 - term_fun2)^2
#	end
#	
#	return total
#end

