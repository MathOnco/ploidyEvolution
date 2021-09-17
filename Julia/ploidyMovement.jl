using Pkg;Pkg.activate(".");Pkg.instantiate();

using DifferentialEquations,ComponentArrays
using LinearAlgebra
using DelimitedFiles,CSV,DataFrames
using Parameters
using Distributed
using Printf



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
function ploidy_wellmixed_model(du,u,pars,t)

	# Grab the parameters
	@unpack (focal_birthRate,nChrom,minChrom,maxChrom,stepChrom,misRate,deathRate) = pars

	# Iterate over each compartment 
	for i in CartesianIndices(u)

		#=
		1) calculateParents(focal, minChrom, maxChrom, stepChrom)
		2) updatePopulation (a struct?)
			sum(inflow(focal,parent,u[parent]) for parent in parentList)
			- outflow(focal,u[focal])
			a. inflow needs to calculate birth and flow rate
			b. outflow needs to calculate birth, death and flow rate
		=#

		focalCN =  collect(i[k] for k in 1:nChrom);

		parentCNList = calculateParents(focalCN, minChrom,maxChrom,stepChrom)

		inflow = 0.0

		# Inflow from the parents
		for parentCN in parentCNList

			# Get birth rate from precomputed
			birthRate = focal_birthRate[Int.(parentCN)...]

			# Get flow rate from parentCN -> focalCN
			flowRate = q(parentCN,focalCN,misRate,nChrom)

			# Will need to have this be the index in the future
			inflow += birthRate*flowRate*u[Int.(parentCN)...] # *(1 - 1/avg(CN))

		end

		# Grab the focal cells birth rate
		birthRate = focal_birthRate[i]

		# Add it to the inflow to the focal compartment
		inflow += birthRate*(1.0 - 2*misRate)*u[i]

		# Death of the focal compartment
		outflow = deathRate*u[i]

		# Update the RHS
		du[i] = inflow - outflow
		
	end

end

function ploidy_spatial_model(du,u,pars,t)

	# Grab the parameters
	@unpack (focal_birthRate,nChrom,chromArray,misRate,deathRate,Γ,Γₑ,ϕ,Ξ,χ,δ,Np,k,E_vessel,
	domain_Dict,boundary_Dict,max_cell_cycle_duration,max_birthRate,debugging) = pars

	# FIX ME:: Assume that all chroms have same lower, upper and step size
	minChrom,stepsize,maxChrom = Int.([chromArray[1] |> f for f in (first,step,length)])

	# Convert chrom array to integer
	intChromArray = [Int(first(x)):Int(step(x)):Int(last(x)) for x in chromArray]
	s, E = u.s, u.E

	coordsList = [1:num_points for num_points in Np]
	dp = (Np .- 1).^(-1)
	sum_dp_invsq = sum(dp.^(-2))

	# Iterate over each coordinate
	for coord in Iterators.product((1:length(coords) for coords in coordsList)...)

		# Get whether we are interior (0), at boundary of vessel (1) or inside vessel (2)
		domain_info = get(domain_Dict,coord,0)

		if domain_info == 2 # Inside blood vessel we make no changes
			du.s[intChromArray...,coord...] .= 0.0
			du.E[coord...] = 0.0
			continue
		end

		# Gets indicies used to calculate spatial impact (assume periodic)
		cardinal_point_info = get_cardinal_gridpoints([elem for elem in coord],Np)

		# Gets the points to be used for spatial grid resolution
		E_cardinal = [E[info...] for info in cardinal_point_info]
		E_minus,E_plus = E_cardinal[end-1:-2:1],E_cardinal[end:-2:1]
		E_focal = E[coord...]

		if domain_info == 1	# boundary of vessel
			normal_vector = boundary_Dict[coord] 	# Grab normal vector
			for i = 1 : length(dp)
				if normal_vector[i] < 0 # The rightward point is in the vessel so we replace it with the B.C.
					E_plus[i] = E_minus[i] - 2.0*dp[i]*k*(E_focal - E_vessel)
				elseif normal_vector[i] > 0 # The leftward point is in the vessel
					E_minus[i] =  E_plus[i] - 2.0*dp[i]*k*(E_focal - E_vessel)
				end
			end
		end

		# Iterate over each compartment
		for focal in Iterators.product((1:length(chr) for chr in chromArray)...)

			# focal cells copy number state
			focalCN = collect(chromArray[k][focal[k]] for k in 1:nChrom)

			# All possible chromosome states
			parentCNList = calculateParents(focalCN,minChrom,maxChrom,stepsize)

			# Initialize inflow to 0
			inflow = 0.0

			# copy state to local variables for brevity
			s_focal = s[focal...,coord...]
			s_cardinal = [s[focal...,info...] for info in cardinal_point_info]
			s_minus,s_plus = s_cardinal[end-1:-2:1],s_cardinal[end:-2:1]

			# Checking the boundary
			if domain_info == 1	# boundary of vessel
				normal_vector = boundary_Dict[coord] 	# Grab normal vector
				for i = 1 : length(dp)
					if normal_vector[i] < 0 # The rightward point is in the vessel so we replace it with the B.C.
						s_plus[i] = s_minus[i] + χ/Γ*s_focal*(chemotaxis_form(E_plus[i],Ξ) - 
						chemotaxis_form(E_minus[i],Ξ))
					elseif normal_vector[i] > 0 # The leftward point is in the vessel
						s_minus[i] = s_plus[i] - χ/Γ*s_focal*(chemotaxis_form(E_plus[i],Ξ) - 
						chemotaxis_form(E_minus[i],Ξ))
					end
				end
			end

			energy_constraint = energy_constrained_birth(E_focal,ϕ)
			
			if max_birthRate*energy_constraint > log(2)/max_cell_cycle_duration
				# Inflow from the parents
				for parentCN in parentCNList

					s_parent = s[Int.(parentCN)...,coord...]

					# Get max birth rate
					birthRate = focal_birthRate[Int.(parentCN)...]

					if s_parent == 0.0 || birthRate == 0.0
						continue
					end

					# If birth rate is below a threshold determined by amount spent in G1 (default is 100 hours)
					# then we assume that birth does not occur...
					energy_modified_birth_rate = birthRate*energy_constraint
					if energy_modified_birth_rate > log(2)/max_cell_cycle_duration
						# Get flow rate from parentCN -> focalCN
						flowRate = q(parentCN,focalCN,misRate,nChrom)
						inflow += energy_modified_birth_rate*flowRate*s_parent
					else
						# println("At energy level $(E_focal), birth rate of $(parentCN) is too low...")
					end

				end


				if s_focal > 0
					# Grab the focal cells (max?) birth rate
					birthRate = focal_birthRate[Int.(focalCN)...]
					energy_modified_birth_rate = birthRate*energy_constraint

					if energy_modified_birth_rate > log(2)/max_cell_cycle_duration
						# Add it to the inflow to the focal compartment
						inflow += energy_modified_birth_rate*(1.0 - 2.0*misRate)*s_focal
					else
						# println("At energy level $(E_focal), birth rate of $(focalCN) is too low...")
					end
				end

			end

			# Death of the focal compartment
			outflow = deathRate*s_focal

			# u_xx -> (u[i+1] + u[i-1] - 2u[i])/dx^2 
			# random migration
			diffusion = Γ*(sum((s_plus[i] + s_minus[i])/dp[i]^2 for i in 1:length(dp)) 
			- 2.0*sum_dp_invsq*s_focal)

			# chemotaxis
			chemotaxis = χ*(
				(
				sum(
					(chemotaxis_form(E_plus[i],Ξ) + chemotaxis_form(E_minus[i],Ξ))/dp[i]^2
					for i in 1 : length(dp)
				)
				- 2.0*sum_dp_invsq*chemotaxis_form(E_focal,Ξ))*s_focal + 
				sum(
					(s_plus[i] - s_minus[i])*(chemotaxis_form(E_plus[i],Ξ) - chemotaxis_form(E_minus[i],Ξ))/dp[i]^2
					for i in 1 : length(dp)))/4.0

			# Update the RHS dependent on where we are in the domain
			du.s[focal...,coord...] = inflow - outflow + diffusion - chemotaxis
		end

		# Diffusion of energy molecule through the tissue
		diffusion = Γₑ*(sum((E_plus[i] + E_minus[i])/dp[i]^2 for i in 1:length(dp)) - 2.0*sum_dp_invsq*E_focal)

		# Consumption by the cells at the grid point designated by coord.
		consumption = δ*E_focal*sum(s[intChromArray...,coord...])

		# update the energy compartment
		du.E[coord...] = diffusion - consumption
		
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

function simulate_wellmixed(params,X::AbstractArray,Y::AbstractVector,
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
	compartmentMinimum,
	progress_check,
	interpolation_order) = params

	# Polyharmonic interpolator
	interp = PolyharmonicInterpolation.PolyharmonicInterpolator(X,Y,interpolation_order)

	# Number of chromosomes (inferred from interp)
	nChrom = interp.dim

	# The allowable states
	chromArray = [minChrom:stepsize:maxChrom for _ in 1:nChrom]

	# Total number of states
	nComp = prod(length(arr) for arr in chromArray)

	# Time interval to run simulation and initial condition
	tspan = (0.0,finalDay)
	u0 = zeros(Int(maxChrom)*ones(Int,nChrom)...)
	u0[startingPopCN...] = 1.0

	##compute birth rates associated with each copy number state
	focal_birthRate = zeros(Int(maxChrom)*ones(Int,nChrom)...);
	for focal in Iterators.product((1:length(j) for j in chromArray)...)
		focalCN = collect(chromArray[k][focal[k]] for k in 1:nChrom);
		focal_birthRate[focal...] = max(PolyharmonicInterpolation.polyharmonicSpline(interp,focalCN')[1],0.0);
	end

	# If we are replating, we do so when the population is million-fold in size
	cb1 = nothing
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
		cb1 = ContinuousCallback(condition,affect!,nothing)

	end

	##### This is printed to terminal to tell you what time point you're at #####
	cb2 = nothing
	if progress_check
		print_time(integrator) = @printf("t = %.2f\n",integrator.t)
		cb2 = PeriodicCallback(print_time, 0.05,save_positions=(false,false))
	end
	cbset = CallbackSet(cb1,cb2)

	if debugging > 0
		println("Beginning simulation...")
	end

	# run simulation
	odePars = (focal_birthRate=focal_birthRate,
				nChrom=nChrom,
				minChrom = chromArray[1].offset, 
				maxChrom = chromArray[1].len, 
				stepChrom=Int(chromArray[1].step),
				misRate=misRate,
				deathRate=deathRate,
				)
	prob = ODEProblem(ploidy_wellmixed_model,u0,tspan,odePars)
	sol = solve(prob,Tsit5(),maxiters=1e5,abstol=1e-8,reltol=1e-5,saveat=1,callback=cbset)

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

	return soln,cnArray,sol.t

end

function simulate_spatial(params,X::AbstractArray,Y::AbstractVector,
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
	Γ,Γₑ,ϕ,Ξ,χ,δ,Np,
	compartmentMinimum,
	progress_check,
	max_cell_cycle_duration,
	interpolation_order) = params

	# convert max cell cycle to days (the time scale used)
	max_cell_cycle_duration/= 24.0

	k,E_vessel,saveat = 0.5,1.0,0.1
	#= Build domain matrix:

	0 : In domain
	1 : Boundary of blood vessel
	2 : Interior of blood vessel

	=#
	df = CSV.read("13496_2_Slides_and_Data_xy_test.txt",DataFrame)
	maxX,maxY=maximum(df[!,"Centroid X µm"]),maximum(df[!,"Centroid Y µm"])
	minX,minY=minimum(df[!,"Centroid X µm"]),minimum(df[!,"Centroid Y µm"])
	df[!,"Centroid X µm"] .= (df[!,"Centroid X µm"] .- minX)./(maxX .- minX)
	df[!,"Centroid Y µm"] .= (df[!,"Centroid Y µm"] .- minY)./(maxY .- minY)

	x = range(0,1,length=Np[1])
	y = range(0,1,length=Np[2])
	dp = (Np .- 1).^(-1)

	df[!,"Centroid X µm"] .= x[searchsortedfirst.(Ref(x),df[!,"Centroid X µm"])]
	df[!,"Centroid Y µm"] .= y[searchsortedfirst.(Ref(y),df[!,"Centroid Y µm"])]
	df = combine(groupby(df,["Centroid X µm","Centroid Y µm"]),last)

	domain_Dict = Dict{Tuple{Int,Int},Int}() # 1 = boundary blood vessel, 2 = inside blood vessel

	for i = 1 : Np[1]
		xdist = (x[i] .- df[!,"Centroid X µm"]).^2
		for j = 1 : Np[2]
			ydist = (y[j] .- df[!,"Centroid Y µm"]).^2
			if minimum(xdist .+ ydist) < sum(dp.^2)/4
				domain_Dict[(i,j)]=1
			end
		end
	end

	# Fill dict with info about the normal vector
	boundary_Dict = Dict{Tuple{Int,Int},Any}()

	# But first we find whether the point is interior
	for k in keys(domain_Dict)
		i,j = k
		i_plus = (i == Np[1]) ? 1 : i+1
		i_minus = (i == 1) ? Np[1] : i-1
		j_plus = (j == Np[2]) ? 1 : j+1
		j_minus = (j == 1) ? Np[2] : j-1

		# Vessel interior point
		if issubset([(i,j_plus),(i,j_minus),(i_plus,j),(i_minus,j)],keys(domain_Dict))
			domain_Dict[k] = 2
		end
	end

	# If no cardinal direction has an interior point we will delete it from the dict
	for (k,v) in domain_Dict

		# We are at an interior point
		if v == 1
			i,j = k[1],k[2]
			i_plus = (i == Np[1]) ? 1 : i+1
			i_minus = (i == 1) ? Np[1] : i-1
			j_plus = (j == Np[2]) ? 1 : j+1
			j_minus = (j == 1) ? Np[2] : j-1

			cardinal_directions = [(i,j_plus),(i,j_minus),(i_plus,j),(i_minus,j)]

			if any(z->z==2,get(domain_Dict,k,0) for k in cardinal_directions)
				# Vessel boundary point
				normal = [0;0]
				normal[1] += get(domain_Dict,(i,j_plus),0) == 2 ? -1 : 0
				normal[1] += get(domain_Dict,(i,j_minus),0) == 2 ? 1 : 0
				normal[2] += get(domain_Dict,(i_plus,j),0) == 2 ? -1 : 0
				normal[2] += get(domain_Dict,(i_minus,j),0) == 2 ? 1 : 0
				# normal/=norm(normal)
				boundary_Dict[k] = normal
			else
				delete!(domain_Dict,k)
			end
		end
	end



	# Initialize energy to be zero outside the interior of the blood vessels
	E0 = zeros(Np...)
	[E0[k...]=E_vessel for (k,v) in domain_Dict if v==2]

	# Polyharmonic interpolator
	interp = PolyharmonicInterpolation.PolyharmonicInterpolator(X,Y,interpolation_order)

	# Number of chromosomes (inferred from interp)
	nChrom = interp.dim

	# The allowable states
	chromArray = [minChrom:stepsize:maxChrom for i in 1:nChrom]

	# Total number of copy number states
	nComp = prod(length(arr) for arr in chromArray)

	# Time interval to run simulation and initial condition
	tspan = (0.0,finalDay)
	s0 = zeros(Int(maxChrom)*ones(Int,nChrom)...,Np...)
	s0[startingPopCN...,:,:] .= rand(Np...)*startPop 	# CN is randomly placed on the grid
	[s0[startingPopCN...,k...] = 0 for k in keys(domain_Dict)]

	u0 = ComponentArray(s=s0,E=E0)

	# If we are replating, we do so when the population is million-fold in size
	cb1 = nothing
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
		cb1 = ContinuousCallback(condition,affect!,nothing)

	end

	##### This is printed to terminal to tell you what time point you're at ######
	##### Add the smallest nonzero compartment size (and index) to the progress check
	cb2 = nothing
	if progress_check
		print_time(integrator) = @printf("t = %.2f\n",integrator.t)
		cb2 = PeriodicCallback(print_time, 0.1,save_positions=(false,false))
	end
	cbset = CallbackSet(cb1,cb2)

	if debugging > 0
		println("Beginning simulation...")
	end

	##compute birth rates associated with each copy number state
	focal_birthRate = zeros(Int(maxChrom)*ones(Int,nChrom)...);
	for focal in Iterators.product((1:length(j) for j in chromArray)...)
		focalCN = collect(chromArray[k][focal[k]] for k in 1:nChrom);
		focal_birthRate[focal...] = max(PolyharmonicInterpolation.polyharmonicSpline(interp,focalCN')[1],0.0);
		println("CN = $focalCN, birth rate = $(focal_birthRate[focal...])")
	end

	max_birthRate = maximum(focal_birthRate)

	# run simulation
	# isoutofdomain = (u,p,t)->any(x->x<0,u)
	odePars = (focal_birthRate=focal_birthRate,
				nChrom=nChrom,
				chromArray=chromArray,
				misRate=misRate,
				deathRate=deathRate,
				Γ=Γ,Γₑ=Γₑ,ϕ=ϕ,Ξ=Ξ,χ=χ,δ=δ,Np=Np,k=k,E_vessel=E_vessel,
				domain_Dict=domain_Dict,
				boundary_Dict=boundary_Dict,
				max_cell_cycle_duration=max_cell_cycle_duration,
				max_birthRate=max_birthRate,
				debugging=debugging)
	prob = ODEProblem(ploidy_spatial_model,u0,tspan,odePars)
	sol = solve(prob,VCABM(),abstol=1e-8,reltol=1e-5,saveat=saveat,callback=cbset)

	# Maybe lsoda() instead of VCABM() ????
	#=

	t = range(tspan..., length=1000) [u.s for u in sol(t)]
	[u.s for u in sol.u]
	=#

	if debugging > 0
		println("Simulation complete. Collecting results...")
	end

	# Grab copy number state array -> CN array
	cnArray = Array{Int}(undef,nComp,nChrom)
	for (i,focal) in enumerate(
		Iterators.product((1:length(j) for j in chromArray)...))
			cnArray[i,:] = collect(chromArray[k][focal[k]] for k in 1:nChrom)
	end

	return sol, cnArray

end

#############################################
#											#
#											#
#					Main		 			#
#											#
#											#
#############################################

function runPloidyMovement(params,X::AbstractArray,Y::AbstractVector,
	startingPopCN::AbstractArray,spatial::Bool)

	if spatial
		println("Running the spatial model...")
		return simulate_spatial(params,X,Y,startingPopCN)
	else
		println("Running the well-mixed model...")
		return simulate_wellmixed(params,X,Y,startingPopCN)
	end

end