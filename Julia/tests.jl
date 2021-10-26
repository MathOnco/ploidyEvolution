include("stochastic.jl")
# define some variables that could plausibly come from deterministic model and
# instantiate stochastic compartment model.Allow stochastic model to run independently 
#of spatial model for 600 timesteps & write output.

function test1()
    #Random.seed!(2)
    Np = [10,10]
    Npde = 2 ## number of active states in the PDE model
    nChrom= 3
    ϕ= 0.05
    dt=0.1
    Γ=.5
    χ=0.1
    Ξ=0.1
    misrate=0.2
    deathRate=0.15
    maxChrom=3


    stochInitPars = (Np=Np,ϕ=ϕ, nChrom=nChrom, maxChrom=maxChrom,dt=dt,
    misrate=misrate,deathRate=deathRate,Γ=Γ, χ=χ,Ξ=Ξ)
    s=setup_stochastic(stochInitPars)

    # assume there will be 2 arrays in scope when the stochastic model is used:
    E0 = zeros(Np...) # energy array
    for i in 1:Np[1]
        for j in 1:Np[2]
            E0[i,j] = 0.5*(i+j)/Np[1]
        end
    end

    

    s0 = zeros(Npde,Np...) # deterministic state array
    s0[1,:,:] .= rand(Np...)*5000
    s0[2,:,:] .= rand(Np...)

    sIndex=([1,1,1],[2,2,2])

    # pre-allocated birth rates for deterministic component will be available
    birthRates=[get_birthrate(sIndex[1]), get_birthrate(sIndex[2])] 
    stochInfo=(sIndex=sIndex,birthRates=birthRates)
        
    for i in 1:600
        s, u =stochastic_step(s,E0,s0,stochInfo)
        println(sum(u))
        if i > 1 
            #can see how the stochastic model evolves by itself by wiping
            # the deterministic array
            s0 = zeros(Npde,Np...) 
        end
        if i%5==0
            df = mapreduce(vcat, values(s.popDict)) do f
                summarise_clone(f)
            end
            saveid = string(i,pad=3)
            CSV.write(saveid*".csv", df)
        end     
    end
end

test1()