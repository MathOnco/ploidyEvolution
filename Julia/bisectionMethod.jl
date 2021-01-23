module Bisection

using LinearAlgebra
import Base.@kwdef

@kwdef struct Results{T<:Real}
	best::T				# The best x that solves the optimization problem
	val::T				# The best value returned (e.g. val = f(best))
	iter::Int			# Number of iterations used in the method
	f_calls::Int
	info::String		# Converged?
end

@kwdef mutable struct Options
	iter::Int
	f_tol::Float64
	tol::Float64
end

function Options(;iter::Int=1000,f_tol::Float64 = 1e-3,tol::Float64 = 1e-3)

	Options(iter,f_tol,tol)

end

function bisectionSearch(a::Real,b::Real,f::Function,options::Options=Options())

	# Promote the endpoints to the same type
	a,b = promote(a,b)

	# Default (not converged)
	converged = false

	# get initial values for high and low 
	fa,fb,f_calls = f(a),f(b),2

	# If they have the same sign then this method will not work
	if sign(fa) == sign(fb)
		@warn("sgn(f(a)) = sgn(f(b)) is not allowed for bisection.")
		return a,fa,b,fb
	end

	# Introduce midpoint and f(midpoint) to be returned to the results struct
	local iterConverged,c, fc

	# Begin bisection
	for iter = 1 : options.iter

		# Compute midpoint
		c = (a + b)/2.0

		# Evaluate f
		fc, f_calls = f(c), f_calls+1

		# Check which sign f(c) has and replace the endpoint with the same sign
		if sign(fa) == sign(fc)
			fa = fc
			a = c
		else
			fb = fc
			b = c
		end

		# Converged
		if abs(fc) < options.f_tol || (b - a) < options.tol
			iterConverged = iter
			converged = true
			break
		end

	end

	if converged
		info = string("Converged in ", iterConverged, " iterations.")
		return Results(c,fc,iterConverged,f_calls,info)
	else
		info = string("Did not converge in ", options.iter, " iterations.")
		return Results(c,fc,options.iter,f_calls,info)
	end

	

end


end



