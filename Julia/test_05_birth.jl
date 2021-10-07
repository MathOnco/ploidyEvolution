using DifferentialEquations,ComponentArrays

## test whether birth and death work properly in the models.
## to achieve this we will generate a 
include("stochastic.jl")
include("ploidyMovement.jl")

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
    starting_copy_number=[2,2,2,2,2]
    maxChrom=5
    misRate=0.0
    deathRate=0.1
    Γ=100.0
    Γₑ=0.0
    ϕ=0.02
    Ξ=0.02
    χ=0.0
    δ=0.01
    Np=[21,21]
    Lp=[1000.0,1000.0]
    dp = Lp./(Np .- 1)
    k=0.0
    E_vessel=0.0
    max_cell_cycle_duration=4
    debugging=0
    stochasticThreshold=1000
    dt=0.1
    tspan = (0.0,dt)
    finalDay=100
    domain_Dict = Dict{Tuple{Int,Int},Int}()
    boundary_Dict = Dict{Tuple{Int,Int},Any}()
    
    CN_birthrate_df = extract_XY("birthLandscapeBrainCancer.txt","GrowthRate_brain_cancer_CLs.txt")
    copy_number,birth_rate = Matrix(CN_birthrate_df[!,r"chr"]),Matrix(CN_birthrate_df[!,r"birth"])
    chromosome_header = names(CN_birthrate_df[!,r"chr"])
    birth_rate = dropdims(birth_rate,dims=2)
    interp = PolyharmonicInterpolation.PolyharmonicInterpolator(copy_number,birth_rate,2)
    nChrom = interp.dim
    
    s0 = zeros(1,Np...) 
    s0[1,:,:] .= 1 	
    E0 = zeros(Np...)
    E0.+=10
    u0 = ComponentArray(s=s0,E=E0)
    
    ##compute birth rates associated with each copy number state
    sIndex=[starting_copy_number]
    birthRates = [max(PolyharmonicInterpolation.polyharmonicSpline(interp,starting_copy_number')[1],0.0)]
    M=transition_matrix(sIndex,birthRates)
    max_birthRate = maximum(birthRates)
    
    stochInitPars = (Lp=Lp,Np=Np,ϕ=ϕ, nChrom=nChrom, maxChrom=maxChrom,dt=dt,
    misrate=misRate,deathRate=deathRate,Γ=Γ, χ=χ,Ξ=Ξ,interp=interp,
    max_cell_cycle_duration=max_cell_cycle_duration, domain_Dict=domain_Dict,
    stochasticThreshold=stochasticThreshold)
    s_s=setup_stochastic(stochInitPars)
    u_s=zeros(Np...)
    
    odePars = (nChrom=nChrom,
    misRate=misRate,
    deathRate=deathRate,
    Γ=Γ,Γₑ=Γₑ,ϕ=ϕ,Ξ=Ξ,χ=χ,δ=δ,Np=Np,dp=dp,k=k,E_vessel=E_vessel,
    domain_Dict=domain_Dict,
    boundary_Dict=boundary_Dict,
    max_cell_cycle_duration=max_cell_cycle_duration,
    max_birthRate=max_birthRate,
    debugging=debugging)
    out=hcat("time","stoch","pde")
    open("test_output/05_birth/cells.csv","w") do io
        writedlm(io,out,',')
    end
    for i in 1:ceil(Int64,finalDay/dt)
        t=i*dt
        s_new, sIndex, birthRates, relegation_index, u_s, append_state = run_hybrid_step(i,odePars,u0,tspan,s_s,sIndex,birthRates,M,u_s)
        # works even if both events happen in same timestep, since new clones are appended thus not affecting the index to be removed
        max_birthRate = maximum(birthRates)
        println(length(sIndex))
        if append_state 
            #println(size(s_new))
            u0= ComponentArray(s=s_new,E=u0.E)
            M=transition_matrix(sIndex,birthRates)
            compartment_sizes = map(x->sum(u0.s[x,:,:]),1:length(sIndex))
            #println(compartment_sizes)
            #println(birthRates)
            #println(sIndex)
        end
        if  relegation_index[1]< stochasticThreshold && sum(u0.s)>0
            ## handle special case when there is still a dummy population in the PDE model
            if sum(sIndex[1])<1
                ## do not add deleted dummy population to the stochastic model
            else
                ## extract info about deleted population
                ui = s_new[relegation_index[2],:,:] 
                br = birthRates[relegation_index[2]]
                cn = sIndex[relegation_index[2]]
                s_s=addClone(s_s, ui, cn, br)
            end
            ## remove state from PDE model
            s_new = s_new[1:size(s_new,1) .!= relegation_index[2],:,:]
            birthRates = birthRates[1:size(birthRates,1) .!= relegation_index[2]]
            sIndex = sIndex[1:size(sIndex,1) .!= relegation_index[2]]
            u0= ComponentArray(s=s_new,E=u0.E)
            if length(sIndex)<1
                println("pde compartment deactivated...")
                ## if no populations above threshold still need to run pde model to handle 
                ## diffusion. PDE function expects a cell population, so we feed it an array of zeros so it will still run.
                u0= ComponentArray(s= zeros(1,Np...) ,E=u0.E)
                sIndex=[repeat([0],length(starting_copy_number))]
                birthRates=0
            end
            compartment_sizes = map(x->sum(u0.s[x,:,:]),1:length(sIndex))
            #println(compartment_sizes)
        end
        
        Nstoch=sum(u_s)
        Npde=sum(u0.s)
    
        out=hcat(t,Nstoch,Npde)
        sE0 = sum(u.E)
        println("time: $t, PDE: $Npde, Stochastic: $Nstoch, Energy $sE0")
        open("test_output/05_birth/cells.csv","a") do io
            writedlm(io,out,',')
        end
    end
    return 0
end

main()