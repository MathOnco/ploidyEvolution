using TOML, CSV, DataFrames, ArgParse, JSON

# Used to search through command line for particular arguments
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input","-i"
            help = "Specifies an inputfile to replace default parameters"
        "--verbosity", "-v"
            help = "Prints information to command line for debugging"
            action = :store_true
		"--spatial", "-s"
			help = "Use the spatial ploidy model"
			action = :store_true
		"--cnFile", "-c"
			help = "Contains numeric array of copy numbers and header of which chromosomes"
			default = "birthLandscapeBrainCancer.txt"
		"--birthRateFile", "-b"
			help = "Contains numeric vector of birth rates"
			default = "GrowthRate_brain_cancer_CLs.txt"
		"--outputfile","-o"
			help = "Specifies the filename for saving data from the simulation"
			default = "output"
		"--u0"
			help = "Array that contains CN of initial condition"
			default = "[1,1,1,1,1]"
    end
	
    return parse_args(s)
end

# Struct containing input data needed for simulation with default values
struct WellMixedInput

	debugging::Int					# prints info from ploidyMovement
	stepsize::Real					# discretization of chromosome
	minChrom::Real					# minimum chromosome allowed
	maxChrom::Real					# maximum chromosome allowed
	deathRate::Float64				# universal death rate
	misRate::Float64 				# universal missegregation rate
	finalDay::Real					# end of simulation
	replating::Bool					# Whether we replate the cells
	startPop::Real					# Starting population size 
	maxPop::Real					# Max population before replating
	compartmentMinimum::Bool		# Sets sizes < 1 to 0 if true
	progress_check::Bool			# Prints the current time of the simulation
	interpolation_order::Int		# Allows the user to set the order of polyharmonic spline

end

struct SpatialInput

	debugging::Int					# prints info from ploidyMovement
	stepsize::Number				# discretization of chromosome
	minChrom::Number
	maxChrom::Number				# maximum chromosome allowed
	deathRate::Number				# universal death rate
	misRate::Float64				# universal missegregation rate
	Γ::Number						# Random cell motion
	Γₑ::Number						# Nutrient diffusion
	ϕ::Number						# Nutrient level needed for half-maximal growth rate
	Ξ::Number						# Energy needed for half maximal directional motion
	χ::Number						# Directed motion magnitude
	δ::Number						# Energy consumption
	Np::Vector{Int}					# Grid size
	finalDay::Real					# end of simulation
	replating::Bool					# Whether we replate the cells
	startPop::Real					# Starting population size 
	maxPop::Real					# Max population before replating
	compartmentMinimum::Bool		# Sets sizes < 1 to 0 if true
	max_cell_cycle_duration::Number	# Length of cell cycle (assumed to be in hours for now)
	progress_check::Bool			# Prints the current time of the simulation
	interpolation_order::Int		# Allows the user to set the order of polyharmonic spline

end

function WellMixedInput()

	# Get the parameters for the struct.
	debugging= 0
	stepsize= 1.0 
	minChrom=1.0
	maxChrom=5.0
	deathRate=0.1
	misRate=0.15
	finalDay=30.0
	replating=false
	startPop=1e3
	maxPop=1e6
	compartmentMinimum=false
	progress_check=false
	interpolation_order=2

	WellMixedInput(
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
		progress_check,
		interpolation_order
		)

end

function WellMixedInput(inputFile::String)

	# Dictionary that contains info from the toml (INPUT) file
	data = TOML.tryparsefile(inputFile)

	# Get the parameters for the struct.
	debugging=get(data,"debugging",0)
	stepsize=get(data,"stepsize",1.0)
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
	interpolation_order = get(data,"interpolation_order",2)
	

	WellMixedInput(
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
		progress_check,
		interpolation_order
		)

end

function SpatialInput()

	# Get the parameters for the struct.
	debugging= 0
	stepsize= 1.0 
	minChrom=1.0
	maxChrom=5.0
	deathRate=0.1
	misRate=0.15
	finalDay=30.0
	replating=false
	startPop=1e3
	maxPop=1e6
	compartmentMinimum=false
	Γ = 0.05
	Γₑ = 0.5
	ϕ = 0.5
	Ξ = 0.5
	χ = 2.0
	δ = 0.005
	Np = [21,21]
	progress_check=false
	interpolation_order=2
	max_cell_cycle_duration = 100.0

	SpatialInput(
		debugging,
		stepsize,
		minChrom,
		maxChrom,
		deathRate,
		misRate,
		Γ,
		Γₑ,
		ϕ,
		Ξ,
		χ,
		δ,
		Np,
		finalDay,
		replating,
		startPop,
		maxPop,
		compartmentMinimum,
		max_cell_cycle_duration,
		progress_check,
		interpolation_order
		)

end

function SpatialInput(inputFile::String)

	# Dictionary that contains info from the toml (INPUT) file
	data = TOML.tryparsefile(inputFile)

	# Get the parameters for the struct.
	debugging=get(data,"debugging",0)
	stepsize=get(data,"stepsize",1.0)
	minChrom = 1.0
	maxChrom=get(data,"maxChrom",5.0)
	deathRate=get(data,"deathRate",0.1)
	misRate=get(data,"misRate",0.15)
	finalDay=get(data,"finalDay",30.0)
	replating=get(data,"replating",false)
	startPop=get(data,"startPop",1e3)
	maxPop=get(data,"maxPop",1e6)
	compartmentMinimum=get(data,"compartmentMinimum",false)
	Γ = get(data,"Γ",0.05)
	Γₑ = get(data,"Γₑ",0.05)
	ϕ = get(data,"ϕ",1.0)
	Ξ = get(data,"Ξ",1.0)
	χ = get(data,"χ",2.0)
	δ = get(data,"δ",0.005)
	Np = get(data,"Np",[21,21])
	progress_check = get(data,"progress_check",false)
	interpolation_order = get(data,"interpolation_order",2)
	max_cell_cycle_duration = get(data,"max_cell_cycle_duration",100.0)

	SpatialInput(
		debugging,
		stepsize,
		minChrom,
		maxChrom,
		deathRate,
		misRate,
		Γ,
		Γₑ,
		ϕ,
		Ξ,
		χ,
		δ,
		Np,
		finalDay,
		replating,
		startPop,
		maxPop,
		compartmentMinimum,
		max_cell_cycle_duration,
		progress_check,
		interpolation_order
		)

end

function initialize()

	# Grab command line args
	parsed_args = parse_commandline()

	# Prints info if set to true
	verbosity = parsed_args["verbosity"]

	spatial = parsed_args["spatial"]

	# Grab file containing copy number variation
	CNmatFilename = parsed_args["cnFile"]

	# Grab file containing birth rates
	birthFilename = parsed_args["birthRateFile"]

	# Grab output file name
	outputFile = parsed_args["outputfile"]

	# Grab initial condition (it is a string so we parse and convert)
	u0 = Int.(JSON.parse(parsed_args["u0"]))

	# Default values to be fed into ploidyMovement
	if spatial
		data = SpatialInput()
		# Replaces default parameters if input file is given
		if !isnothing(parsed_args["input"])
			data = SpatialInput(parsed_args["input"])
		end
	else
		data = WellMixedInput()
		# Replaces default parameters if input file is given
		if !isnothing(parsed_args["input"])
			data = WellMixedInput(parsed_args["input"])
		end
	end

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

	return data, CNmatFilename, birthFilename, u0, outputFile, verbosity, spatial
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

end
	

function main()

	# Grab input parameters and whether to print info to terminal
	data, CNmatFilename, birthFilename, u0, outputFile, verbosity, spatial = initialize()

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
	# results, time = runPloidyMovement(data,copy_number,birth_rate,u0,spatial)
	
	if spatial

		sol, cnArray = runPloidyMovement(data,copy_number,birth_rate,u0,spatial)

		# Write results to multiple files by time points
		for (index,t) in enumerate(sol.t)
			open(outputFile*"_states_time_"*replace(string(t),"." => "_")*".csv","w") do io
				for (i,cnstate) in enumerate(eachrow(cnArray))
					z = sol[index].s[cnstate...,:,:]
					cnstate = "#"*join(string(cnArray[i,:]...),":")
					writedlm(io,[cnstate])
					writedlm(io,z,',')
				end
			end
			open(outputFile*"_Energy_time_"*replace(string(t),"." => "_")*".csv","w") do io
				writedlm(io,sol[index].E,',')
			end
		end

		# Save discretized domain
		writedlm(outputFile*"_x.csv",range(0,1,length=data.Np[1]))
		writedlm(outputFile*"_y.csv",range(0,1,length=data.Np[2]))

	else

		sol, cnArray,time = runPloidyMovement(data,copy_number,birth_rate,u0,spatial)

		# convert to array for output
		results = hcat(cnArray,sol)
		outputHeader = permutedims(vcat(chromosome_header,time))

		# concatnate the header with the results array
		output = vcat(outputHeader,results)
	
		# save to file
		writedlm( outputFile*".csv",  output, ',')

	end

end

@time main()
