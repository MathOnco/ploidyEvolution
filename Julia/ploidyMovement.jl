using Base: Number
using Pkg;Pkg.activate(".");Pkg.instantiate();

using DifferentialEquations,ComponentArrays
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
	(interp,nChrom,chromArray,misRate,deathRate,migrationRate,Np,debugging) = pars

	##### FIXME ######
	intChromArray = [Int(first(x)):Int(step(x)):Int(last(x)) for x in chromArray]
	s, E = u.s, u.E
	ϕ = 1.0
	Ξ = 1.0
	χ = 3.0
	consumption_rate = 0.001
	energy_diffusion_rate = 0.05

	# @show maximum.(chromArray)

	coordsList = [1:num_points for num_points in Np]
	dp = (Np .- 1).^(-1)
	sum_dp_invsq = sum(dp.^(-2))

	# Iterate over each compartment
	for coord in Iterators.product((1:length(coords) for coords in coordsList)...)

		# Points used to calculate spatial impact (assume periodic)
		cardinal_point_info = get_cardinal_gridpoints([elem for elem in coord],Np)

		E_cardinal = [E[info...] for info in cardinal_point_info]
		E_minus,E_plus = E_cardinal[1:2:end],E_cardinal[2:2:end]
		Etmp = E[coord...]

		for focal in Iterators.product((1:length(chr) for chr in chromArray)...)

	#=
	1) calculateParents(focal, minChrom, maxChrom, stepChrom)
	2) updatePopulation (a struct?)
		sum(inflow(focal,parent,u[parent]) for parent in parentList)
		- outflow(focal,u[focal])
		a. inflow needs to calculate birth and flow rate
		b. outflow needs to calculate birth, death and flow rate
	=#

			focalCN = collect(chromArray[k][focal[k]] for k in 1:nChrom)

			parentCNList = calculateParents(focalCN, 1,2,1)		##### FIXME #####

			inflow = 0.0

			# copy state to local variables for brevity
			stmp = s[focal...,coord...]
			s_cardinal = [s[focal...,info...] for info in cardinal_point_info]
			s_minus,s_plus = s_cardinal[1:2:end],s_cardinal[2:2:end]
			

			# Inflow from the parents
			for parentCN in parentCNList

				# Get (max?) birth rate
				birthRate = max(PolyharmonicInterpolation.polyharmonicSpline(interp,parentCN')[1],0.0)

				# Get flow rate from parentCN -> focalCN
				flowRate = q(parentCN,focalCN,misRate,nChrom)

				# Will need to have this be the index in the future
				inflow += birthRate*energy_constrained_birth(Etmp,ϕ)*flowRate*s[Int.(parentCN)...,coord...]

			end

			# Grab the focal cells (max?) birth rate
			birthRate = max(PolyharmonicInterpolation.polyharmonicSpline(interp,focalCN')[1],0.0)

			# Add it to the inflow to the focal compartment
			inflow += birthRate*energy_constrained_birth(Etmp,ϕ)*(1.0 - 2*misRate)*stmp

			# Death of the focal compartment
			outflow = deathRate*stmp

			# u_xx -> (u[i+1] + u[i-1] - 2u[i])/dx^2 
			# random migration
			diffusion = migrationRate*(sum((s_plus[i] + s_minus[i])/dp[i]^2 for i in 1:length(dp)) 
			- 2.0*sum_dp_invsq*stmp)

			# chemotaxis
			chemotaxis = χ*((sum((chemotaxis_form(E_plus[i],Ξ) + chemotaxis_form(E_plus[i],Ξ))/dp[i]^2 for i in 1 : length(dp))
			- 2.0*sum_dp_invsq*chemotaxis_form(Etmp,Ξ))*stmp + 
			sum( (s_plus[i] - s_minus[i])*(chemotaxis_form(E_plus[i],Ξ) - chemotaxis_form(E_minus[i],Ξ))/dp[i]^2
			for i in 1 : length(dp)))/4.0

			# Information if debugging
			if debugging > 5 && (inflow != 0.0 || outflow != 0.0)
				@show t,focalCN,inflow,outflow, s[focal...]
			end

			# Update the RHS
			du.s[focal...,coord...] = inflow - outflow + diffusion - chemotaxis
		end

		diffusion = energy_diffusion_rate*(sum((E_plus[i] + E_minus[i])/dp[i]^2 for i in 1:length(dp)) 
			- 2.0*sum_dp_invsq*Etmp)

		consumption = consumption_rate*Etmp*sum(s[intChromArray...,coord...])

		# update the energy
		du.E[coord...] = diffusion - consumption
		
	end



	if debugging > 3
		@show t,sum(s)
	end

end

function get_cardinal_gridpoints(coordinate::Vector{Int},Np::Vector{Int})

	@assert length(coordinate) == length(Np) "Dimensions must match"
	dim = length(coordinate)

	[ (tmp = copy(coordinate); tmp[idx] = mod(v-1,Np[idx])+1; tmp) 
		for idx in 1:dim for v in coordinate[idx]-1:coordinate[idx]+1 
			if v != coordinate[idx]]
end

function calculateParents(offspring::Vector{T}, minChrom::Int,
	maxChrom::Int,stepChrom::Int) where T<:Real


	[ (x = copy(offspring); x[idx] = v; x)
		for idx in 1 : length(offspring)
			for v in max(offspring[idx]-1,minChrom):stepChrom:min(offspring[idx]+1,maxChrom) 
				if v != offspring[idx] ]

end

function energy_constrained_birth(energy::Number,ϕ::Number)

	return energy/(ϕ + energy)

end

function chemotaxis_form(energy::Number,Ξ::Number)

	return log(Ξ + energy)

end

function energy_consumption(energy::Number,δ::Number,u::Vector{T}) where T <: Number

	return δ*energy*sum(u)

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

	Np = [21,21]
	migrationRate = 0.005
	E0 = ones(Np...)
	E0[1:10,:] .= 0.0

	# Polyharmonic interpolator
	interp = PolyharmonicInterpolation.PolyharmonicInterpolator(X,Y)

	# Number of chromosomes (inferred from interp)
	nChrom = interp.dim

	# minimum and maximum size of copy number
	minX,maxX = minimum(X,dims=1),maximum(X,dims=1)

	# The allowable states
	chromArray = [minChrom:stepsize:maxChrom for i in 1:nChrom]

	# Total number of states
	nComp = prod(length(arr) for arr in chromArray)*prod(Np)

	# Time interval to run simulation and initial condition
	tspan = (0.0,finalDay)
	s0 = zeros(Int(maxChrom)*ones(Int,nChrom)...,Np...)
	s0[startingPopCN...,:,:] .= rand(Np...)						# CN is uniformly placed on the grid

	u0 = ComponentArray(s=s0,E=E0)

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
	odePars = (interp,nChrom,chromArray,misRate,deathRate,migrationRate,Np,debugging)
	prob = ODEProblem(ploidyModel,u0,tspan,odePars)
	sol = solve(prob,Tsit5(),maxiters=1e5,abstol=1e-8,reltol=1e-5,saveat=1,callback=callback)

	#=

	t = range(tspan..., length=1000) [u.s for u in sol(t)]
	[u.s for u in sol.u]
	=#

	if debugging > 0
		println("Simulation complete. Collecting results...")
	end

	return sol

	# # convert to array for output
	# soln = reshape(Array(sol),nComp,length(sol.t))

	# # Grab compartment -> CN array
	# cnArray = Array{Float64}(undef,nComp,nChrom)
	# for (i,focal) in enumerate(
	# 	Iterators.product((1:length(j) for j in chromArray)...))
	# 		cnArray[i,:] = collect(chromArray[k][focal[k]] for k in 1:nChrom)
	# end

	# # Concatnate the chrom indices with soln output
	# results = hcat(cnArray,soln)

	# return results, sol.t

end