using Pkg;Pkg.activate(".");Pkg.instantiate();

using TOML, DataFrames, CSV, ArgParse, Distributed, Dates

@everywhere import Base.@kwdef

# Used to search through command line for particular arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--paramFile","-p"
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
		"--initfile", "-i"
			help = "file that contains simulation information"
    end

    return parse_args(s)
end

# Struct containing input data needed for simulation with default values
@everywhere @kwdef struct Input

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

@everywhere function Input(paramFile::String)

	# Dictionary that contains info from the toml (INPUT) file
	params = TOML.tryparsefile(paramFile)

	# Get the parameters for the struct.
	debugging=get(params,"debugging",0)
	stepsize=get(params,"stepsize",1.0)
	minChrom=get(params,"minChrom",1.0)
	maxChrom=get(params,"maxChrom",5.0)
	deathRate=get(params,"deathRate",0.1)
	misRate=get(params,"misRate",0.15)
	finalDay=get(params,"finalDay",30.0)
	replating=get(params,"replating",false)
	startPop=get(params,"startPop",1e3)
	maxPop=get(params,"maxPop",1e6)
	compartmentMinimum=get(params,"compartmentMinimum",false)

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
	params = Input()

	# Grab command line args
	parsed_args = parse_commandline()

	# Replaces default parameters if input file is given
	if !isnothing(parsed_args["paramFile"]) && isfile(parsed_args["paramFile"])
		params = Input(parsed_args["paramFile"])
	end

	# Prints info if set to true
	verbosity = parsed_args["verbosity"]

	# Grab file containing copy number variation
	CNmatFilename = parsed_args["cnFile"]

	# Grab file containing birth rates
	birthFilename = parsed_args["birthRateFile"]

	# Grab output file name
	outputFile = parsed_args["outputfile"]

	# Grab initial condition file
	initfile = parsed_args["initfile"]

	# Error checking
	if !isfile(CNmatFilename)
		error("$CNmatFilename file cannot be found")
	end
	if !isfile(birthFilename)
		error("$birthFilename file cannot be found")
	end

	# Check to see if input file is given
    if verbosity
		println("parsed_args:")
		for (arg,val) in parsed_args
			println("  $arg => $val")
		end
	end

	return params, CNmatFilename, birthFilename, initfile, outputFile, verbosity
end

function main()

	# Grab input parameters and whether to print info to terminal
	params, CNmatFilename, birthFilename, initfile, outputFile, verbosity = initialize()

	date = today() # Get today's date to add to outputfile string

	# Printing stuff to terminal
	if verbosity
		println("Input:")
		for name in fieldnames(typeof(params)) 
			println("  $name => $(getfield(params,name))")
		end
	end

	# Read the birth rate file to be used for the interpolation
	X = readdlm(CNmatFilename, '\t')
	birthRates = readdlm(birthFilename, '\t')

	# FIXME!!!!!
	header, cn = X[1,:],Float64.(X[2:end,3:end])

	# Error if the elements in X or birthRates are not all subtypes of real
	if any((eltype(cn),eltype(birthRates)) .>: Real)
		error("cn and birthRates must contain only numeric values")
	end

	# FIXME!!!!!!! At the moment we assume ploidy is not included
	cn .+= 2.0 # adding ploidy by hand.

	# readdlm gives a 2D array, we turn it into a vector (1d array) here.
	birthRates = dropdims(birthRates,dims=2)

	# Load initial condition from file if available:
	if !isfile(initfile)
		error("paralleldriver must have an initfile")
	else
		df = CSV.read(initfile,DataFrame)
	end

	# Grab the initial conditions and convert to integer for now
	u0mat = Int.(round.(convert(Matrix,df[:,[name for name in names(df)
	if occursin(r"Chr*",name)]])))

	# check if we update either of the two rates
	changeDeathRate = "deathRate" in names(df) ? true : false
	changeMisRate = "missegregationRate" in names(df) ? true : false

	if verbosity
		println("Running simulation in parallel with $(nworkers()) workers")
	end

	# Parallelization of the initial condition file (u0mat)
	pmap(enumerate(eachrow(u0mat))) do (i,initCN)

		# Define the outputfile string
		outputFile = "output_$(date)_run-$i.csv"

		# Check to see if we replace the death or misRate
		if changeDeathRate
			deathRate = df[i,"deathRate"]
		else
			deathRate = -1.0
		end
		if changeMisRate
			misRate = df[i,"missegregationRate"]
		else
			misRate = -1.0
		end

		# run the simulation
		results, time = runPloidyMovement(params,cn,birthRates,initCN,deathRate,
		misRate)

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

	end
	
end

main()
