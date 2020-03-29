module Utils

export Resolution
export get_sinewave
export get_evolved
export get_evolved_slice
export commandlineparser
export commandlineparser_int
export check_time

using Burgers

function get_sinewave(x)
	# returns a sine wave with a period if 1 on [0, 1]
	return 0.35 .+ 0.25*sin.(2*pi*x)
end

struct Resolution{T, I}
	# contains information about the resolution of the 2D (x, t) grid
	# that will be used by the FunArr2D instances
	Ttot::T
	Nt::I
	Ncells::I
	Nghost::I
end

function check_time(res::Resolution{T}, t0::T) where {T}
	# makes sure that the given time is valid for the resolution
	return 0.0 <= t0 < res.Ttot
end

function get_evolved_slice(func::Function, res::Resolution{T, I}, t::T, sc::SlopeCalc)::FunArr{T, T} where {T, I<:Int}
	# evaluates 'func' on [0, 1] to initialize u(x, 0), and evolves in time,
	# then returns the function u(x) evaluated at a time t
	it = Int64(res.Nt*floor(t/res.Ttot))
	
	init_positions = FunArr{Float64, Float64}(func, res.Ncells, res.Nghost)
	u2fun          = evolve(init_positions, res.Nt, res.Ttot, sc, it)
	ufun_t         = get_timeslice(u2fun, t, res.Ncells, res.Nghost)
	
	return ufun_t
end

function get_evolved_slice(func::Function, res::Resolution{T, I}, it::I, sc::SlopeCalc)::FunArr{T, T} where {T, I<:Int}
	# evaluates 'func' on [0, 1] to initialize u(x, 0), and evolves in time,
	# then returns the function u(x) evaluated at the timeslice it
	init_positions = FunArr{Float64, Float64}(func, res.Ncells, res.Nghost)
	u2fun          = evolve(init_positions, res.Nt, res.Ttot, sc, it)
	uarr_t         = u2fun.u2dat[:,it]
	ufun_t         = FunArr{T, T}(uarr_t, res.Ncells, 0)
	
	return ufun_t
end

function get_evolved(func::Function, res::Resolution{T, I}, sc::SlopeCalc)::FunArr2D{T, T} where {T, I<:Int}
	# evaluates 'func' on [0, 1] to initialize u(x, 0), and evolves in time,
	# then returns the function u(x) evaluated at a time t
	init_positions = FunArr{Float64, Float64}(func, res.Ncells, res.Nghost)
	u2fun          = evolve(init_positions, res.Nt, res.Ttot, sc, res.Nt)
	
	return u2fun
end

function commandlineparser(check_int::Bool)::Tuple{Number, SlopeCalc, AbstractString}
	# parses the command line arguments
	if length(ARGS) != 3
		if check_int == true
			error("This program requires 3 arguments: (1) index for the evolved time, (2) the type of slope calculation to perform, and (3) the name of the file where the image will be saved.")
		else
			error("This program requires 3 arguments: (1) the evolved time, (2) the type of slope calculation to perform, and (3) the name of the file where the image will be saved.")
		end
	end
	
	t_str, sc_str, fname = ARGS
	
	# parse the slope calculation type
	sc =
	if sc_str == "minmod"
		minmod::SlopeCalc
	elseif sc_str == "linear"
		linear::SlopeCalc
	elseif sc_str == "constant"
		constant::SlopeCalc
	else
		error("the second argument must be one of: minmod linear constant")
	end
	
	# parse the time
	if check_int == true
		it = tryparse(Int64, t_str)
		if typeof(it) == Nothing
			error("the first argument must be an Int64 value.")
		end
		return it, sc, fname
	else
		t0 = tryparse(Float64, t_str)
		if typeof(t0) == Nothing
			error("the first argument must be a Float64 value.")
		end
		return t0, sc, fname
	end
end

end #module
