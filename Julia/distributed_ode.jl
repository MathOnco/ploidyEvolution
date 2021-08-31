using Pkg;Pkg.activate(".");Pkg.instantiate();
using TOML, DelimitedFiles, ArgParse, JSON,  DifferentialEquations, BenchmarkTools, Distributed, Parameters
import Base.@kwdef



# Load the file that contains the Polyharmonic interpolator function
include("polyHarmonicInterp.jl")

# Used to search through command line for particular arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input","-i"
            help = "Specifies an inputfile to replace default parameters"
        "--verbosity", "-v"
            help = "Prints information to command line for debugging"
            action = :store_true
		"--cnFile", "-c"
			help = "Contains numeric array of copy numbers and header of which chromosomes"
			default = "birthLandscapeBrainCancer.txt"
		"--birthRateFile", "-b"
			help = "Contains numeric vector of birth rates"
			default = "myY.txt"
		"--outputfile","-o"
			help = "Specifies the filename for saving data from the simulation"
			default = "output_0.csv"
		"--u0"
			help = "Array that contains CN of initial condition"
			default = "[1,1,1,1,1]"
    end
	
    return parse_args(s)
end

# Struct containing input data needed for simulation with default values
@kwdef struct Input

	debugging::Int = 0				# prints info from ploidyMovement
	stepsize::Real = 1.0				# discretization of chromosome
	minChrom::Real = 1.0				# minimum chromosome allowed
	maxChrom::Real = 5.0				# maximum chromosome allowed
	deathRate::Float64 = 0.1			# universal death rate
	misRate::Float64 = 0.15				# universal missegregation rate
	finalDay::Real = 30.0				# end of simulation
	replating::Bool = false				# Whether we replate the cells
	startPop::Real = 1e3		# Starting population size 
	maxPop::Real = 1e6				# Max population before replating
	compartmentMinimum::Bool = false	# Sets sizes < 1 to 0 if true

end

function Input(inputFile::String)

	# Dictionary that contains info from the toml (INPUT) file
	data = TOML.tryparsefile(inputFile)

	# Get the parameters for the struct.
	debugging=get(data,"debugging",0)
	stepsize=get(data,"stepsize",1.0)
	minChrom=get(data,"minChrom",1.0)
	maxChrom=get(data,"maxChrom",5.0)
	deathRate=get(data,"deathRate",0.1)
	misRate=get(data,"misRate",0.15)
	finalDay=get(data,"finalDay",30.0)
	replating=get(data,"replating",false)
	startPop=get(data,"startPop",1e3)
	maxPop=get(data,"maxPop",1e6)
	compartmentMinimum=get(data,"compartmentMinimum",false)

	Input(
		debugging,
		stepsize,
		minChrom,
		maxChrom,
		deathRate,
		misRate,
		finalDay,
		replating,
		startPop,
		maxPop,
		compartmentMinimum
		)

end

function initialize()

	# Default values to be fed into ploidyMovement
	data = Input()

	# Grab command line args
	parsed_args = parse_commandline()

	# Replaces default parameters if input file is given
	if !isnothing(parsed_args["input"])
		data = Input(parsed_args["input"])
	end

	# Prints info if set to true
	verbosity = parsed_args["verbosity"]

	# Grab file containing copy number variation
	CNmatFilename = parsed_args["cnFile"]

	# Grab file containing birth rates
	birthFilename = parsed_args["birthRateFile"]

	# Grab output file name
	outputFile = parsed_args["outputfile"]

	# Grab initial condition (it is a string so we parse and convert)
	u0 = Int.(JSON.parse(parsed_args["u0"]))

	# Error checking
	if !isfile(CNmatFilename)
		error("$CNmatFilename file cannot be found")
	end
	if !isfile(birthFilename)
		error("$birthFilename file cannot be found")
	end

	# Avoid overwriting an old .csv file if using default
	count = 1
	while isfile(outputFile)
		if verbosity
			println("outpileFile = $outputFile exists, changing to default names.")
		end
		outputFile = "output_$count.csv"
		count += 1
	end

	# Check to see if input file is given
    if verbosity
		println("parsed_args:")
		for (arg,val) in parsed_args
			println("  $arg => $val")
		end
	end

	return data, CNmatFilename, birthFilename, u0, outputFile, verbosity
end

function ploidyModel_para(du,u,pars,t)

	#du_shared= SharedArray(du)
	# Grab the parameters
	(interp,nChrom,chromArray,misRate,deathRate,debugging,birthRates) = pars

	# Iterate over each compartment
    @sync @distributed for i in CartesianIndices(u)
        inflow=0.0
        outflow=0.0

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
		birthRate = birthRates[i]#max(PolyharmonicInterpolation.polyharmonicSpline(interp,focalCN')[1],0.0)

		# Add it to the inflow to the focal compartment
		inflow += birthRate*(1.0 - 2*misRate)*u[i]

		# Death of the focal compartment
		outflow = deathRate*u[i]

		# Update the RHS
		du[i] = inflow - outflow
		
	end

    #du = sdata(du_shared)

end




## main execution here

# Grab input parameters and whether to print info to terminal
data, CNmatFilename, birthFilename, u0, outputFile, verbosity = initialize()
println("initialized")
# Printing stuff to terminal
if verbosity
    println("Input:")
    for name in fieldnames(typeof(data)) 
        println("  $name => $(getfield(data,name))")
    end
end

# Read the birth rate file to be used for the interpolation
X = readdlm(CNmatFilename, '\t')
Y = readdlm(birthFilename, '\t')

header, cn = X[1,:],Float64.(X[2:end,3:end])

# Error if the elements in X or Y are not all subtypes of real
if any((eltype(cn),eltype(Y)) .>: Real)
    error("cn and Y must contain only numeric values")
end

# readdlm gives a 2D array, we turn it into a vector (1d array) here.
Y = dropdims(Y,dims=2)
println("entering runPloidyMovement")

## rename args of what was ploidyMovement
params=data
X=cn
startingPopCN=u0

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
odePars = (interp,nChrom,chromArray,misRate,deathRate,debugging,birthRates)
prob = ODEProblem(ploidyModel_para,u0,tspan,odePars)

addprocs(4)

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
@everywhere using SharedArrays
@everywhere using LinearAlgebra

println("starting ODE solver")
sol = solve(prob,Tsit5(),maxiters=1e5,abstol=1e-8,reltol=1e-5,saveat=1,callback=callback,save_everystep=false)

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

time = sol.t

# create header for the solution output
outputHeader = permutedims(
    vcat(
        ["Chr$chr" for chr in header if typeof(chr) == Int],
        time
        )
        )

# concatnate the header with the results array
output = vcat(outputHeader,results)

# save to file
writedlm( outputFile,  output, ',')

print("done")
