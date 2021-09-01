#using Pkg;Pkg.activate(".");Pkg.instantiate();

using DifferentialEquations, BenchmarkTools, Distributed

using DelimitedFiles
using Parameters

# Load the file that contains the Polyharmonic interpolator function
include("polyHarmonicInterp.jl")

println("init ploidyMovement")

@everywhere using SharedArrays
@everywhere using LinearAlgebra


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
	(maxChrom,interp,nChrom,chromArray,misRate,deathRate,debugging,birthRates) = pars

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
			birthRate = birthRates[Int.(parentCN)...]# max(PolyharmonicInterpolation.polyharmonicSpline(interp,parentCN')[1],0.0)

			# Get flow rate from parentCN -> focalCN
			flowRate = q(parentCN,focalCN,misRate,nChrom)

			# Will need to have this be the index in the future
			inflow += birthRate*flowRate*u[Int.(parentCN)...] # *(1 - 1/avg(CN))

		end

		# Grab the focal cells birth rate
		birthRate = birthRates[focal...]#max(PolyharmonicInterpolation.polyharmonicSpline(interp,focalCN')[1],0.0)

		# Add it to the inflow to the focal compartment
		inflow += birthRate*(1.0 - 2*misRate)*u[focal...]

		# Death of the focal compartment
		outflow = deathRate*u[focal...]

		# Information if debugging
		if debugging > 5 && (inflow != 0.0 || outflow != 0.0)
			@show t,focalCN,inflow,outflow, u[focal...]
		end

		# Update the RHS
		du[focal...] = inflow - outflow
		
	end

	if debugging > 3
		@show t,sum(u)
	end

end

function ploidyModel_para(u,pars,t)

	# Grab the parameters
	(maxChrom, interp,nChrom,chromArray,misRate,deathRate,debugging,birthRates) = pars
	du = SharedArray(zeros(Int(maxChrom)*ones(Int,nChrom)...))

	# Iterate over each compartment
    @sync @distributed for i in CartesianIndices(u)
		inflow = 0.0
		#=
		1) calculateParents(focal, minChrom, maxChrom, stepChrom)
		2) updatePopulation (a struct?)
			sum(inflow(focal,parent,u[parent]) for parent in parentList)
			- outflow(focal,u[focal])
			a. inflow needs to calculate birth and flow rate
			b. outflow needs to calculate birth, death and flow rate
		=#

		focalCN =  collect(i[k] for k in 1:nChrom) 
		parentCNList = calculateParents(focalCN, 1,5,1)

		parent_birthrates = map(xx -> birthRates[Int.(xx)...], parentCNList);
		parent_flowrates = map(xx -> q(xx,focalCN,misRate,nChrom), parentCNList) ;
		parent_u = map(xx -> u[Int.(xx)...], parentCNList);

		inflow += sum(parent_birthrates .* parent_flowrates .* parent_u)

		# Inflow from the parents
	#	for parentCN in parentCNList

			# Get birth rate
		#	birthRate = birthRates[Int.(parentCN)...]# max(PolyharmonicInterpolation.polyharmonicSpline(interp,parentCN')[1],0.0)

			# Get flow rate from parentCN -> focalCN
		#	flowRate = q(parentCN,focalCN,misRate,nChrom)

			# Will need to have this be the index in the future
		#	inflow += birthRate*flowRate*u[Int.(parentCN)...] # *(1 - 1/avg(CN))

		#end

		# Grab the focal cells birth rate
		#birthRate = birthRates[i]#max(PolyharmonicInterpolation.polyharmonicSpline(interp,focalCN')[1],0.0)

		# Add it to the inflow to the focal compartment
		inflow += birthRates[i] .* (1.0 .- 2 .* misRate) .* u[i]

		# Death of the focal compartment
		#outflow = deathRate*u[i]

		# Update the RHS
		#inflow[i] = inflow_i
		du[i]=inflow .- deathRate .* u[i]
		
	end

	return(du)

end


@everywhere function calculateParents(offspring::Vector{T}, minChrom::Int,
	maxChrom::Int,stepChrom::Int) where T<:Real


	[ (x = copy(offspring); x[idx] = v; x)
		for idx in 1 : length(offspring)
			for v in max(offspring[idx]-1,minChrom):stepChrom:min(offspring[idx]+1,maxChrom) 
				if v != offspring[idx] ]

end

@everywhere function q(parent::Vector{T},offspring::Vector{T},misRate::Float64,nChrom::Int) where T <: Real

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
	startingPopCN::AbstractArray)

	#addprocs(4)
	#println("currently running processes: "+nprocs())

	# Grab the parameters from the struct
	@unpack (debugging,
	stepsize,
	minChrom,
	maxChrom,
	deathRate,
	misRate,
	finalDay,
	replating,
	startPop,
	maxPop,
	compartmentMinimum) = params

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
	du = copy(u0)
	u0[startingPopCN...] = 1.0

	birthRates = zeros(Int(maxChrom)*ones(Int,nChrom)...)

	for (i,focal) in enumerate(
		Iterators.product((1:length(j) for j in chromArray)...))
		focalCN = collect(chromArray[k][focal[k]] for k in 1:nChrom)
		birthRate = max(PolyharmonicInterpolation.polyharmonicSpline(interp,focalCN')[1],0.0)
		birthRates[focal...] = birthRate
	end
	# If we are replating, we do so when the population is million-fold in size
	callback = nothing
	if replating
		# Replate if above maxPop
		condition = (u,t,integrator) -> sum(u) - maxPop	

		# This replates so that u sums to 1 (as in the initial condition)
		affect!(integrator) = begin
			if debugging > 1
				println("Max population reached at t = $(integrator.t)...replating.")
			end
			integrator.u /= maxPop
		end

		# Create callback
		callback = ContinuousCallback(condition,affect!,nothing)

	end

	# If we artificially set populations below threshold to 0
	# if compartmentMinimum
	# 	condition = (u,t,integrator)

	if debugging > 0
		println("Beginning simulation...")
	end

	#include("dist_functions.jl")
	#addprocs(4)

	# run simulation
	odePars = (maxChrom,interp,nChrom,chromArray,misRate,deathRate,debugging,birthRates)
	prob = ODEProblem(ploidyModel,u0,tspan,odePars)
	println("starting ODE solver")
	sol = solve(prob,Tsit5(),maxiters=1e5,abstol=1e-8,reltol=1e-5,saveat=1,callback=callback)

	if debugging > 0
		println("Simulation complete. Collecting results...")
	end

	# convert to array for output
	soln = reshape(Array(sol),nComp,length(sol.t))

	# Grab compartment -> CN array
	cnArray = Array{Float64}(undef,nComp,nChrom)
	for (i,focal) in enumerate(
		Iterators.product((1:length(j) for j in chromArray)...))
			cnArray[i,:] = collect(chromArray[k][focal[k]] for k in 1:nChrom)
	end

	# Concatnate the chrom indices with soln output
	results = hcat(cnArray,soln)

	# # create header for the solution output
	# header = permutedims(vcat(["chr$i" for i in 1 : nChrom],sol.t))

	# # concatnate the header with the results array
	# output = vcat(header,results)

	# # save to file
	# writedlm( outputFile,  output, ',')

	return results, sol.t

end