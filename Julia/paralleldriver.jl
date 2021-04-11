using Pkg;Pkg.activate(".");Pkg.instantiate();

using TOML, DelimitedFiles, ArgParse, JSON, Distributed, Dates

@everywhere import Base.@kwdef

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
                "--u0file"
                    help = "file that contains initial condition"
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

@everywhere function Input(inputFile::String)

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

	# Grab initial condition file
	u0file = parsed_args["u0file"]

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

	return data, CNmatFilename, birthFilename, u0file, outputFile, verbosity
end

function main()

	# Grab input parameters and whether to print info to terminal
	data, CNmatFilename, birthFilename, u0file, outputFile, verbosity = initialize()

	date = today() # Get today's date to add to outputfile string

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

	# FIXME!!!!!
	header, cn = X[1,:],Float64.(X[2:end,3:end])

	# Error if the elements in X or Y are not all subtypes of real
	if any((eltype(cn),eltype(Y)) .>: Real)
		error("cn and Y must contain only numeric values")
	end

	# FIXME!!!!!!! At the moment we assume ploidy is not included
	cn .+= 2.0 # adding ploidy by hand.

	# readdlm gives a 2D array, we turn it into a vector (1d array) here.
	Y = dropdims(Y,dims=2)

	# Load initial condition from file if available:
	if !isfile(u0file)
		if verbosity
			println("Initial condition file NOT found. using default")
		end
		u0 = [2,2,2,2,2]
	else
		u0mat = readdlm(u0file,'\t')
	end

	# Parallelization of the initial condition file (u0mat)
	pmap(enumerate(eachrow(u0mat))) do (i,u0)
		ploidy,cnv = u0[1],u0[2:end]
		initCN = Int.(round.(ploidy .+ cnv))
		outputFile = "output_$(date)_run-$i.csv"
		results, time = runPloidyMovement(data,cn,Y,initCN)

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
