using Pkg;Pkg.activate(".");Pkg.instantiate();

using DifferentialEquations
using LinearAlgebra
using DelimitedFiles
using Parameters

# Load the file that contains the Polyharmonic interpolator function
include("polyHarmonicInterp.jl")


"""

	plodyModel()

	This function is used in DifferentialEquations and returns the RHS of the
	differential equation needed to solve ploidy evolution of the discrete 
	system.

	INPUTS
		du:		A vector of the RHS of y' = du
		u:		A vector of the solution variable
		p:		A list of parameters that we pass
		t:		The current time point
	OUTPUTS
		du:		A vector of the RHS of y' = du


"""
function ploidyModel(du,u,pars,t)

	# Grab the parameters
	(interp,nChrom,chromArray,misRate,deathRate,debugging) = pars

	# Iterate over each compartment
	for (i,focal) in enumerate(
		Iterators.product((1:length(j) for j in chromArray)...))

		#=
		1) calculateParents(focal, minChrom, maxChrom, stepChrom)
		2) updatePopulation (a struct?)
			sum(inflow(focal,parent,u[parent]) for parent in parentList)
			- outflow(focal,u[focal])
			a. inflow needs to calculate birth and flow rate
			b. outflow needs to calculate birth, death and flow rate
		=#

		focalCN = collect(chromArray[k][focal[k]] for k in 1:nChrom)

		parentCNList = calculateParents(focalCN, 1,5,1)

		inflow = 0.0

		# Inflow from the parents
		for parentCN in parentCNList

			# Get birth rate
			birthRate = max(PolyharmonicInterpolation.polyharmonicSpline(interp,parentCN')[1],0.0)

			# Get flow rate from parentCN -> focalCN
			flowRate = q(parentCN,focalCN,misRate,nChrom)

			# Will need to have this be the index in the future
			inflow += birthRate*flowRate*u[Int.(parentCN)...] # *(1 - 1/avg(CN))

		end

		# Grab the focal cells birth rate
		birthRate = max(PolyharmonicInterpolation.polyharmonicSpline(interp,focalCN')[1],0.0)

		# Add it to the inflow to the focal compartment
		inflow += birthRate*(1.0 - 2*misRate)*u[focal...]

		# Death of the focal compartment
		outflow = deathRate*u[focal...]

		# Information if debugging
		if debugging && (inflow != 0.0 || outflow != 0.0)
			@show t,focalCN,inflow,outflow, u[focal...]
		end

		# Update the RHS
		du[focal...] = inflow - outflow
		
	end

end

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

#############################################
#											#
#											#
#					Main		 			#
#											#
#											#
#############################################

function runPloidyMovement(params,X::AbstractArray,Y::AbstractVector,
	startingPopCN::AbstractArray,outputFile::String)

	# Grab the parameters from the struct
	@unpack (debugging,
	stepsize,
	minChrom,
	maxChrom,
	deathRate,
	misRate,
	finalDay) = params

	# Polyharmonic interpolator
	interp = PolyharmonicInterpolation.PolyharmonicInterpolator(X,Y)

	# Number of chromosomes (inferred from interp)
	nChrom = interp.dim

	# minimum and maximum size of copy number
	minX,maxX = minimum(X,dims=1),maximum(X,dims=1)

	# The allowable states
	chromArray = [minChrom:stepsize:maxChrom for i in 1:nChrom]

	# Total number of states
	nComp = prod(length(arr) for arr in chromArray)

	# Time interval to run simulation and initial condition
	tspan = (0.0,finalDay)
	u0 = zeros(Int(maxChrom)*ones(Int,nChrom)...)
	u0[startingPopCN...] = 1.0

	# run simulation
	odePars = (interp,nChrom,chromArray,misRate,deathRate,debugging)
	prob = ODEProblem(ploidyModel,u0,tspan,odePars)
	sol = solve(prob,Tsit5(),maxiters=1e5,abstol=1e-8,reltol=1e-5,saveat=1)

	# convert to array for output
	soln = reshape(Array(sol),nComp,length(sol.t))

	# Grab compartment -> CN array
	cnArray = Array{Float64}(undef,nComp,nChrom)
	for (i,focal) in enumerate(
		Iterators.product((1:length(j) for j in chromArray)...))
			cnArray[i,:] = collect(chromArray[k][focal[k]] for k in 1:nChrom)
	end

	# Concatnate the solution
	output = hcat(cnArray,soln)

	# save to file
	writedlm( outputFile,  output, ',')

	return

end
