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
			default = "GrowthRate_brain_cancer_CLs.txt"
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
	progress_check::Bool = false

end

function Input(inputFile::String)

	# Dictionary that contains info from the toml (INPUT) file
	data = TOML.tryparsefile(inputFile)

	# Get the parameters for the struct.
	debugging=get(data,"debugging",0)
	stepsize=get(data,"stepsize",1.0)
	# minChrom=get(data,"minChrom",1.0)
	minChrom=1.0
	maxChrom=get(data,"maxChrom",5.0)
	deathRate=get(data,"deathRate",0.1)
	misRate=get(data,"misRate",0.15)
	finalDay=get(data,"finalDay",30.0)
	replating=get(data,"replating",false)
	startPop=get(data,"startPop",1e3)
	maxPop=get(data,"maxPop",1e6)
	compartmentMinimum=get(data,"compartmentMinimum",false)
	progress_check=get(data,"progress_check",false)

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
		compartmentMinimum,
		progress_check
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

	# Check to see if input file is given
    if verbosity
		println("parsed_args:")
		for (arg,val) in parsed_args
			println("  $arg => $val")
		end
	end

	return data, CNmatFilename, birthFilename, u0, outputFile, verbosity
end

function extract_XY(X_filename::String, Y_filename::String)

	# Read the birth rate file to be used for the interpolation
	X_df = CSV.read(X_filename,DataFrame)
	Y_df = CSV.read(Y_filename, DataFrame)

	# We make all strings upper case to ensure consistency
	[X_df[!,name] = uppercase.(x) for (name,x) in zip(names(X_df),eachcol(X_df)) 
					if eltype(x) === String]
	[Y_df[!,name] = uppercase.(y) for (name,y) in zip(names(Y_df),eachcol(Y_df)) 
					if eltype(y) === String]

	# Join the data frames by cell line name
	XY_df = outerjoin(X_df, Y_df,on = intersect(names(X_df), names(Y_df)))

	# Error if the number of elements in X_df or Y_df do not match
	@assert nrow(X_df) == nrow(Y_df) == nrow(XY_df) "Cell line names likely did not match."

	# Once we are sure the sizes are the same, we remove missings
	dropmissing!(XY_df)

	return XY_df

	# return X,Y

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

	# Extract copy number (X) and birth rate (Y) to be used in the interpolation in ploidyMovement.jl
	CN_birthrate_df = extract_XY(CNmatFilename,birthFilename)

	# Get copy number and birth rates
	local copy_number,birth_rate,chromosome_header
	try
		copy_number,birth_rate = Matrix(CN_birthrate_df[!,r"chr"]),Matrix(CN_birthrate_df[!,r"birth"])
		chromosome_header = names(CN_birthrate_df[!,r"chr"])
	catch
		@error("chromosome X should be the column names")
	end

	# Drop dims is required to convert to a vector (from mat, required for polyharmonic)
	birth_rate = dropdims(birth_rate,dims=2)

	# Run ploidy movement
	results, time = runPloidyMovement(data,copy_number,birth_rate,u0)

	outputHeader = permutedims(vcat(chromosome_header,time))

	# create header for the solution output
	# outputHeader = permutedims(
	# 	vcat(
	# 		["Chr$chr" for chr in header if typeof(chr) == Int],
	# 		time
	# 		)
	# 		)

	# concatnate the header with the results array
	output = vcat(outputHeader,results)

	# save to file
	writedlm( outputFile,  output, ',')
	
end

@time main()
