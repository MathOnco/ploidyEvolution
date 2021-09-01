using Pkg;Pkg.activate(".");Pkg.instantiate();

using Distributed, DifferentialEquations, DelimitedFiles, Parameters, CUDA
@everywhere using SharedArrays
@everywhere using LinearAlgebra

# Load the file that contains the Polyharmonic interpolator function
include("polyHarmonicInterp.jl")


"""

	ploidyModel()

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
	(interp,nChrom,chromArray,misRate,deathRate,focal_birthRate,compartments,debugging) = pars

	parent_inflow = SharedArray(zeros(Int(chromArray[1].len)*ones(Int,nChrom)...));

	# Iterate over each compartment
	@sync @distributed for (i,focal) in collect(enumerate(compartments))

		#=
		1) calculateParents(focal, minChrom, maxChrom, stepChrom)
		2) updatePopulation (a struct?)
			sum(inflow(focal,parent,u[parent]) for parent in parentList)
			- outflow(focal,u[focal])
			a. inflow needs to calculate birth and flow rate
			b. outflow needs to calculate birth, death and flow rate
		=#
		
		focalCN = collect(chromArray[k][focal[k]] for k in 1:nChrom);

		parentCNList = calculateParents(focalCN, chromArray[1].offset, chromArray[1].len, Int(chromArray[1].step));

		# Get flow rate from parentCN -> focalCN
		flowRate = map(t -> q(t,focalCN,misRate,nChrom), parentCNList) ;
		# flowRate = cu(flowRate_);

		# Get birth rate from parentCN 
		birthRate = map(t -> max(PolyharmonicInterpolation.polyharmonicSpline(interp,t')[1],0.0), parentCNList);
		# birthRate = cu(birthRate_);

		# Get size of parental compartments
		v = map(t -> u[Int.(t)...], parentCNList);
		# v = cu(v_)

		# inflow from parentCN 
		parent_inflow[focal...] = sum(birthRate.*flowRate.*v)

	end
	#foci= cu(u);

	# Add parental inflow to the inflow to the focal compartment	
	inflow = parent_inflow .+  focal_birthRate.*(1.0 - 2*misRate).*u;

	outflow = deathRate * u;

	du .= Array(inflow .- outflow);

	#print("\nt=",t,": ",maximum(du))
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
	compartments=Iterators.product((1:length(j) for j in chromArray)...);

	# Grab the focal cells birth rate and compartment size
	focal_birthRate = zeros(Int(maxChrom)*ones(Int,nChrom)...);
	for (i,focal) in enumerate(compartments)
		focalCN = collect(chromArray[k][focal[k]] for k in 1:nChrom);
		focal_birthRate[focal...] = max(PolyharmonicInterpolation.polyharmonicSpline(interp,focalCN')[1],0.0);
	end

	# Time interval to run simulation and initial condition
	tspan = (0.0,finalDay)
	# u0 = cu(zeros(Int(maxChrom)*ones(Int,nChrom)...))
	# CUDA.@allowscalar u0[startingPopCN...] = 1.0
	u0 = zeros(Int(maxChrom)*ones(Int,nChrom)...)
	u0[startingPopCN...] = 1.0

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

	# run simulation
	odePars = (interp,nChrom,chromArray,misRate,deathRate,focal_birthRate, compartments, debugging)
	prob = ODEProblem(ploidyModel,u0,tspan,odePars)
	sol = solve(prob,Tsit5(),maxiters=1e5,abstol=1e-8,reltol=1e-5,saveat=1,callback=callback)

	if debugging > 0
		println("Simulation complete. Collecting results...")
	end

	# convert to array for output
	soln = reshape(Array(sol),nComp,length(sol.t))

	# Grab compartment -> CN array
	cnArray = Array{Float32}(undef,nComp,nChrom)
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
