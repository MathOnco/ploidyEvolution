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
		"--outputfile","-o"
			help = "Specifies the filename for saving data from the simulation"
			default = "output"
    end
	
    return parse_args(s)
end

# Convert string keys to symbols (needed for splatting dictionary into struct fields)
string_to_symbol(d::Dict) = Dict(Symbol(k) => v for (k,v) in d)

# Modified get function to get struct from the dict file
function get_structs_from_file(d::Dict,key::String,default::DataType,verbosity::Bool)

	local s
	try
		s = default(string_to_symbol(d[key]))
		if verbosity
			println("Replaced $key with inputfile")
		end
	catch
		s = default()
		if verbosity 
			println("Did not replace $key with inputfile")
		end
	end

	return s

end

struct Options

	debugging::Int
	compartmentMinimum::Bool
	replating::Bool
	maxPop::Number
	progress_check::Bool
	saveat::Number
	interpolation_order::Int

end

struct Files
	copy_number_file::String		# Contains numeric array of copy numbers and header of which chromosomes
	birth_rate_file::String			# Contains numeric vector of birth rates
	blood_vessel_file::String		# Blood vessel coordinate info
end

struct UniversalParameters
	stepsize::Number
	minChrom::Number
	maxChrom::Number
	deathRate::Number
	misRate::Number
	finalDay::Number
	startPop::Number
	starting_copy_number::Vector{Int}
	max_cell_cycle_duration::Number
end

struct SpatialParameters
	Γ::Number						# Random cell motion
	Γₑ::Number						# Nutrient diffusion
	ϕ::Number						# Nutrient level needed for half-maximal growth rate
	Ξ::Number						# Energy needed for half maximal directional motion
	χ::Number						# Directed motion magnitude
	δ::Number						# Energy consumption
	Np::Vector{Int}
	k::Number
	E_vessel::Number
	
end

# Struct containing input data needed for simulation with default values
struct WellMixedInput

	Options
	Files
	UniversalParameters

end

struct SpatialInput

	Options
	Files
	UniversalParameters
	SpatialParameters

end

struct WellMixedInputParameters

	Options
	UniversalParameters

end

struct SpatialInputParameters

	Options
	UniversalParameters
	SpatialParameters
	blood_vessel_file::String

end

function Files(;
	copy_number_file ="birthLandscapeBrainCancer.txt",
	birth_rate_file = "GrowthRate_brain_cancer_CLs.txt",
	blood_vessel_file ="13496_2_Slides_and_Data_xy_test.txt")

	Files(
		copy_number_file,
		birth_rate_file,
		blood_vessel_file
		)
end

Files(d::Dict) = Files(;d...)

function UniversalParameters(;stepsize= 1.0, 
	minChrom=1.0,
	maxChrom=5.0,
	deathRate=0.1,
	misRate=0.15,
	finalDay=30.0,
	startPop=1e3,
	starting_copy_number = [1,1,1,1,1],
	max_cell_cycle_duration = 100)

	UniversalParameters(
					stepsize,
					minChrom,
					maxChrom,
					deathRate,
					misRate,
					finalDay,
					startPop,
					starting_copy_number,
					max_cell_cycle_duration
					)

end

UniversalParameters(d::Dict) = UniversalParameters(;d...)

function Options(;debugging = 3,
	compartmentMinimum = false,
	replating = true,
	maxPop = 1e6,
	progress_check = true,
	saveat = 0.1,
	interpolation_order = 2)

	Options(
			debugging,
			compartmentMinimum,
			replating,
			maxPop,
			progress_check,
			saveat,
			interpolation_order,
			)

end

Options(d::Dict) = Options(;d...)



function SpatialParameters(;
	Γ = 0.1,
	Γₑ = 1.0,
	ϕ = 0.02,
	Ξ = 0.02,
	χ = 0.5,
	δ = 0.05,
	Np = [21,21],
	k = 0.5,
	E_vessel = 1.0)

	SpatialParameters(
					Γ,
					Γₑ,	
					ϕ,
					Ξ,
					χ,
					δ,
					Np,
					k,
					E_vessel
					)



end

SpatialParameters(d::Dict) = SpatialParameters(;d...)

function WellMixedInput(;Options=Options(),
	Files=Files(),
	UniversalParameters=UniversalParameters())

	WellMixedInput(
					Options,
					Files,
					UniversalParameters
					)

end

function WellMixedInput(inputfile::String,verbosity::Bool=false)

	data = TOML.tryparsefile(inputfile)

	WellMixedInput(
					get_structs_from_file(data,"Options",Options,verbosity),
					get_structs_from_file(data,"Files",Files,verbosity),
					get_structs_from_file(data,"UniversalParameters",UniversalParameters,verbosity)
					)


end

function SpatialInput(;Options=Options(),
	Files=Files(),
	UniversalParameters=UniversalParameters(),
	SpatialParameters=SpatialParameters())

	SpatialInput(
					Options,
					Files,
					UniversalParameters,
					SpatialParameters
					)
end

function SpatialInput(inputfile::String,verbosity::Bool=false)

	data = TOML.tryparsefile(inputfile)

	SpatialInput(
					get_structs_from_file(data,"Options",Options,verbosity),
					get_structs_from_file(data,"Files",Files,verbosity),
					get_structs_from_file(data,"UniversalParameters",UniversalParameters,verbosity),
					get_structs_from_file(data,"SpatialParameters",SpatialParameters,verbosity)
					)


end

function initialize()

	# Grab command line args
	parsed_args = parse_commandline()

	# Prints info if set to true
	verbosity = parsed_args["verbosity"]

	spatial = parsed_args["spatial"]

	# Grab output file name
	outputFile = parsed_args["outputfile"]

	# # Grab initial condition (it is a string so we parse and convert)
	# starting_copy_number = Int.(JSON.parse(parsed_args["starting_copy_number"]))

	# Default values to be fed into ploidyMovement
	if spatial
		input = SpatialInput()
		# Replaces default parameters if input file is given
		if !isnothing(parsed_args["input"])
			input = SpatialInput(parsed_args["input"],verbosity)
		end
	else
		input = WellMixedInput()
		# Replaces default parameters if input file is given
		if !isnothing(parsed_args["input"])
			input = WellMixedInput(parsed_args["input"],verbosity)
		end
	end

	# Check to see if input file is given
    if verbosity
		println("parsed_args:")
		for (arg,val) in parsed_args
			println("  $arg => $val")
		end
		println()
	end

	return input, outputFile, verbosity, spatial
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
	input, outputFile, verbosity, spatial = initialize()

	# Printing stuff to terminal
	if verbosity
		for name in fieldnames(typeof(input))
			println("$(name):")
			for var in fieldnames(typeof(getfield(input,name)))
				println("  $var => $(getfield(getfield(input,name),var))")
			end
			println()
		end
	end

	# Error checking
	if !isfile(input.Files.copy_number_file)
		error("$(input.Files.copy_number_file) file cannot be found")
	end
	if !isfile(input.Files.birth_rate_file)
		error("$(input.Files.birth_rate_file) file cannot be found")
	end

	# Extract copy number (X) and birth rate (Y) to be used in the interpolation in ploidyMovement.jl
	CN_birthrate_df = extract_XY(input.Files.copy_number_file,input.Files.birth_rate_file)

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

	# Run ploidy movement and grab params from the input
	### FIX ME: We want to replace input file with dictionaries that we can grab the param chunk without deleting

	if spatial

		params = SpatialInputParameters(input.Options,
										input.UniversalParameters,
										input.SpatialParameters,
										input.Files.blood_vessel_file)

		sol, cnArray = runPloidyMovement(params,copy_number,birth_rate,spatial)

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
		writedlm(outputFile*"_x.csv",range(0,1,length=params.SpatialParameters.Np[1]))
		writedlm(outputFile*"_y.csv",range(0,1,length=params.SpatialParameters.Np[2]))

	else

		params = WellMixedInputParameters(input.Options,input.UniversalParameters)

		sol, cnArray,time = runPloidyMovement(params,copy_number,birth_rate,spatial)

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
