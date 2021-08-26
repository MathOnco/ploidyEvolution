using TOML, CSV, DataFrames, ArgParse, JSON

import Base.@kwdef

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
			default = "GrowthRate brain cancer CLs.txt"
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

<<<<<<< HEAD
=======
	# Avoid overwriting an old .csv file if using default
	count = 1
	while isfile(outputFile)
		if verbosity
			println("outpileFile = $outputFile exists, changing to default names.")
		end
		outputFile = "output_$count.csv"
		count += 1
	end

>>>>>>> a93a09da88bfe7b732a1f273b7217d00b608fa69
	# Check to see if input file is given
    if verbosity
		println("parsed_args:")
		for (arg,val) in parsed_args
			println("  $arg => $val")
		end
	end

	return data, CNmatFilename, birthFilename, u0, outputFile, verbosity
end

function main()

	# Grab input parameters and whether to print info to terminal
	data, CNmatFilename, birthFilename, u0, outputFile, verbosity = initialize()

	# Printing stuff to terminal
	if verbosity
		println("Input:")
		for name in fieldnames(typeof(data)) 
			println("  $name => $(getfield(data,name))")
		end
	end

	# Read the birth rate file to be used for the interpolation
	copy_number_df = CSV.read(CNmatFilename,DataFrame)
	birth_rate_df = CSV.read(birthFilename, DataFrame)

	# We make all strings upper case to ensure consistency
	[copy_number_df[!,name] = uppercase.(x) for (name,x) in zip(names(copy_number_df),eachcol(copy_number_df)) 
					if eltype(x) === String]
	[birth_rate_df[!,name] = uppercase.(y) for (name,y) in zip(names(birth_rate_df),eachcol(birth_rate_df)) 
					if eltype(y) === String]

	# Join the data frames by cell line name
	birth_rate_datapoints = outerjoin(copy_number_df, birth_rate_df, 
								on = intersect(names(copy_number_df), names(birth_rate_df)))

	# Error if the number of elements in copy_number_df or birth_rate_df do not match
	@assert nrow(copy_number_df) == nrow(birth_rate_df) == nrow(birth_rate_datapoints) "Cell line names likely did not match."

	# Once we are sure the sizes are the same, we remove missings
	dropmissing!(birth_rate_datapoints)

	# Get copy number and birth rates
	copy_number,birth_rate = Matrix(birth_rate_datapoints[!,r"chr"]),Matrix(birth_rate_datapoints[!,r"birth"])

	# Drop dims is required to convert to a vector (from mat, required for polyharmonic)
	birth_rate = dropdims(birth_rate,dims=2)

	# Run ploidy movement
	results, time = runPloidyMovement(data,copy_number,birth_rate,u0)

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

main()
