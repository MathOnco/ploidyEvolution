using KrylovKit
using LinearAlgebra
using DelimitedFiles
using Plots

# Load the file that contains the Polyharmonic interpolator function
include("polyHarmonicInterp_v4.jl")
include("bisectionMethod.jl")

function calculateParents(offspring::Vector{T}, minChrom::Int,
	maxChrom::Int,stepChrom::Int) where T<:Real


	[ (x = copy(offspring); x[idx] = v; x)
		for idx in 1 : length(offspring)
			for v in max(offspring[idx]-1,minChrom):stepChrom:
				min(offspring[idx]+1,maxChrom) if v != offspring[idx] ]

end

function q(parent::Vector{T},offspring::Vector{T},misRate::Float64,nChrom::Int) where T <: Real

	# parentCN = collect(linspaces[k][parent[k]] for k in 1:nChrom)
	# offspringCN = collect(linspaces[k][offspring[k]] for k in 1:nChrom)

	parentCN = parent
	offspringCN = offspring

	CNdist = norm(parentCN .- offspringCN,1)

	# No mutation occurred
	if (CNdist == 0.0)
		val = 1.0 - misRate
	else
		# FIXME: Should be a function dependent on parent ploidy
		val = misRate/nChrom
	end

	return val

end

"""

	jacobianMap()
	
	This function defines our Jacobian map which is used by Krylov to calculate
	the maximum eigenvalues of the Jacobian matrix. Of interest is the point at
	which the maximum eigenvalue is 0 (the critical curve) which separates
	extinction from growth.

	INPUTS
		x: 			A vector which is acted on the map of J (J*x -> y)
		mutRate:	The mutation rate of the ploidy types
		birthRate:	The birth rate of the ploidy types
		deathRate:	The death rate of the ploidy types
		matSize:	The number of ploidy types

	OUTPUTS
		y: 			The result of J*x -> y


"""
function jacobianMap(x::Vector{Float64}, misRate::Float64, pars)

	# Initialize the output vector
	y = zeros(length(x))

	(interp,nChrom,chromArray,deathRate,maxY) = pars

	# Linear indexing
	lidx = LinearIndices(zeros(5,5,5,5,5))


	# Iterate over each compartment
	for (i,focal) in enumerate(
		Iterators.product((1:length(j) for j in chromArray)...))

		focalCN = collect(chromArray[k][focal[k]] for k in 1:nChrom)

		parentCNList = calculateParents(focalCN, 1,5,1)

		inflow = 0.0

		# Inflow from other parents
		for parentCN in parentCNList

			# Get birth rate
			birthRate = min(max(PolyharmonicInterpolation.polyharmonicSpline(interp,parentCN')[1],0.0),maxY)

			# Get flow rate from parentCN -> focalCN
			flowRate = q(parentCN,focalCN,misRate,nChrom)

			# Will need to have this be the index in the future
			inflow += birthRate*flowRate*x[lidx[parentCN...]]

		end

		# Grab the focal cells birth rate
		birthRate = min(max(PolyharmonicInterpolation.polyharmonicSpline(interp,focalCN')[1],0.0),maxY)

		inflow += birthRate*(1.0 - 2*misRate)*x[i]

		# Death term
		outflow = deathRate*x[i]

		y[i] = inflow - outflow
		
	end

	return y

end

function getParameterInfo()

	# Read the birth rate file to be used for the interpolation
	X = readdlm("CNVs brain cancer CLs.txt", '\t')
	Y = readdlm("GrowthRate brain cancer CLs.txt", '\t')

	# Prune data to keep only the numeric part
	X = float.(X[2:end,2:end]) .+ 2.0
	Y = float.(Y[:,2])

	# Polyharmonic interpolator
	interp = PolyharmonicInterpolation.PolyharmonicInterpolator(X,Y)

	nChrom = interp.dim

	minX,maxX = minimum(X,dims=1),maximum(X,dims=1)

	stepsize = 1

	chromArray = [1:stepsize:5 for i in 1:nChrom]

	# Total number of compartments
	matSize = prod(length(arr) for arr in chromArray)

	return interp,nChrom,chromArray,matSize,maximum(Y)

end

function getEigenvalues(misRate::Float64,matSize::Int,constParams::Tuple)

	# Find maximum eigenvalues
	evals,evecs,info = (eigsolve(matSize, 1,:LR; verbosity=0,issymmetric=false,
	krylovdim=50,tol=1e-4)
	do u
		# This computes y = J*u (the action of J on the vector u)
		jacobianMap(u,misRate,constParams)
	end)

	# Check if eigenvalue converged
	try
		return real(evals[1])
	catch
		error("Eigenvalue did not converge")
	end

end


# Parameters
deathRateVec = (0.0:0.05:0.7)
interp,nChrom,chromArray,matSize,maxY = getParameterInfo()

# Initialize solution vector
beta = Vector{Float64}()

for deathRate in reverse(deathRateVec)

	# Params for the linear map
	constParams = (interp,nChrom,chromArray,deathRate,maxY)

	# Define the function
	myfun = misRate -> getEigenvalues(misRate,matSize,constParams)

	# Execute bisection search method
	results = Bisection.bisectionSearch(0.0,1.0,myfun)

	# sgn(f(a)) = sgn(f(b)) - bisection does not work. We check to see the sign
	# individually to determine whether any misRate is acceptable or none is.
	if typeof(results) <: Tuple
		a,fa,b,fb = results
		
		# If f(b) < 0 - Any missegregation rate is acceptable
		if fb < 0
			push!(beta,0.0)
		else
			push!(beta,NaN)
		end
	# The signs are unequal and the solver went to work...
	else
		push!(beta,results.best)
	end

	@show deathRate,beta[end]

end



# # Initialize solution vector
# beta = Vector{Float64}()

# for deathRate in reverse(deathRateVec)

# 	result = (optimize(mRate -> lossFunction(mRate,
# 	matSize,(interp,nChrom,chromArray,deathRate)),[0.05], LBFGS()))

# 	@show deathRate, result.minimizer[1]

# 	push!(beta,result.minimizer[1])

# 	initGuess = result.minimizer[1]

# end

# # Beta is a probability and so we remove unphysical solutions by the solver
# # (Alternatively, we can do box-constrained optimization if necessary)
# beta[beta .> 1.0] .= 1.0
# beta[beta .< 0.0] .= 0.0

