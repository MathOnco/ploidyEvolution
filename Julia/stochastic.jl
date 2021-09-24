using DataStructures, Distributions, CSV, DataFrames, Random, Parameters

## TO DO
# make move() failsafe
# input a threshold max number of cells
# increase compatibility with ODE model (e.g. polyharmonicSpline birth rates; use chemotaxis_form() fn's; & etc)
# modify birth_rate() to allow testing with different maxChrom 
mutable struct karyotype
    state 
    viable 
end

function missegregate(k::karyotype, maxChrom,chromosome_selector)
    #select a chromosome randomly
    index = rand(chromosome_selector,1)[1]
    k_new=deepcopy(k)
    #generate daughters according to rule
    k.state[index]=k.state[index]+1
    k_new.state[index]=k_new.state[index]-1
    #check viability
    if k.state[index[1]]<1 || k.state[index[1]]>maxChrom
        k.viable=false
    end
    if k_new.state[index[1]]<1 || k_new.state[index[1]]>maxChrom
        k_new.viable=false
    end
    return k, k_new
end

## handy to set up some sort of forced trajectory of evolution for testing the model
## f(cn) = birthrate: this could easily be swapped with polyharmspline in future
function get_birthrate(cn)
    ## if [1,1,1] is an input state it should go through these states in order:
    if cn == [1,1,1]
        return 0.3
    end
    if cn == [1,1,2]
        return 0.3
    end

    if cn == [1,2,1]
        return 0.3
    end

    if cn == [2,1,1]
        return 0.3
    end

    if cn == [1,1,3]
        return 0.5
    end

    if cn == [1,2,3]
        return 0.6
    end

    if cn == [1,3,3]
        return 0.7
    end

    if cn == [2,3,3]
        return 0.8
    end

    return 0

end

#a cell object
mutable struct cell_allinfo
    x::Float32
    y::Float32
    state
    cell_allinfo(x,y,state) = new(x,y,state) # constructor
end

# a reduced cell object (cells often don't need to know their karyotype)
mutable struct cell
    x::Float32
    y::Float32
    cell(x,y) = new(x,y) # constructor
end

# migrate a cell
function move(c::cell,Np,dnorm, dE0,E0,χ,Ξ)
    ix = ceil(Int64,c.x)
    iy = ceil(Int64,c.y)
    new_x = c.x + rand(dnorm,1)[1] + χ*dE0[1,ix,iy]*log(E0[ix,iy]+Ξ)
    new_y = c.y + rand(dnorm,1)[1] + χ*dE0[2,ix,iy]*log(E0[ix,iy]+Ξ)

    ##implement boundary conditions (NOTE this can be broken if migration longer than Np])
    ## consider % operator to fix overstepping? or just do..while
    if new_x > Np[1] 
        new_x=new_x-Np[1]
    end 
    if new_y > Np[1] 
        new_y=new_y-Np[1]
    end 
    if new_x < 0 
        new_x=new_x+Np[1]
    end 
    if new_y < 0 
        new_y=new_y+Np[1]
    end 
    c.x = new_x
    c.y = new_y
    return c
end

#store all the cells for a given clone
mutable struct clone
    pop
    cn
    birthrate
end

# wrapper function for move()
function migrate(cl::clone,Np,dnorm,dE0,E0,χ,Ξ)
    cl.pop = map(x->move(x,Np,dnorm,dE0,E0,χ,Ξ),cl.pop)
    return(cl)
end

# loop thru all cells per clone, and perform divisions
function replicate(cl::clone, E0, ϕ, dt, uniform01,misrate,maxChrom,chromosome_selector,Np,dnorm,dE0,χ,Ξ)
    # for missegregations we will take a strategy of returning all the new viable karyotypes
    # and then allowing the stochasticCompartment object to deal with them
    missegs = MutableLinkedList{cell_allinfo}()
    for i in length(cl.pop):-1:1
        energy = E0[ceil(Int64,cl.pop[i].x),ceil(Int64,cl.pop[i].y)]
        #need to replace with call to energy function
        probBirth = (cl.birthrate-misrate)*dt*energy/(ϕ + energy)
        if rand(uniform01,1)[1] < probBirth  # division occurs
            if rand(uniform01,1)[1] < misrate #missegregation occurs
                
                k=karyotype(deepcopy(cl.cn),true)
                k1,k2 = missegregate(k,maxChrom,chromosome_selector)
                if k1.viable
                    c1=cell_allinfo(cl.pop[i].x,cl.pop[i].y,k1.state)
                    push!(missegs,c1)
                end
                if k2.viable
                    c2=cell_allinfo(cl.pop[i].x,cl.pop[i].y,k2.state)
                    push!(missegs,c2)
                end
                delete!(cl.pop,i)
            else
                new_c = cell(cl.pop[i].x,cl.pop[i].y)
                new_c=move(new_c,Np,dnorm,dE0,E0,χ,Ξ)
                push!(cl.pop,new_c)
            end
        end
    end
    return cl,missegs
end

# test if all cells per clone need to die, and kill them
function die(cl::clone,dt, deathRate, uniform01)
    # require iterating backwards when deleting or else
    # we will overrun the list

    for i in length(cl.pop):-1:1
        probDeath = deathRate*dt
        if rand(uniform01,1)[1] < probDeath  # division occurs
            delete!(cl.pop,i)
        end
    end
    return(cl)
end

# top level data structure for stochastic compartment, storing state of all clones
# and any other useful variables etc.
mutable struct stochasticCompartment
    popDict # structure to hold all clones (dictionary)
    #store parms used in stochastic model for convenience:
    Np
    ϕ
    χ
    Ξ
    maxChrom
    dt
    misrate
    deathRate
    # store the Distributions used in stochastic model:
    chromosome_selector
    uniform01
    dnorm
end

# perform migration for entire stochastic model
function migrate(s::stochasticCompartment,dE0,E0)
    for k in keys(s.popDict)
        s.popDict[k] = migrate(s.popDict[k],s.Np,s.dnorm,dE0,E0,s.χ,s.Ξ)
    end
    return(s)
end
# perform cell division for entire stochastic model
function replicate(s::stochasticCompartment, E0,dE0,sIndex)
    for k in keys(s.popDict)
        s.popDict[k],missegs = replicate(s.popDict[k],E0, s.ϕ, s.dt, s.uniform01,
        s.misrate,s.maxChrom,s.chromosome_selector,s.Np,s.dnorm,dE0,s.χ,s.Ξ)
        for m in missegs
            if (m.state in sIndex)==false # check if already in det. comp.
                c = cell(m.x,m.y)
                c=move(c,s.Np,s.dnorm,dE0,E0,s.χ,s.Ξ)
                if (m.state in keys(s.popDict)) == false 
                    pop = MutableLinkedList{cell}()
                    push!(pop, c)
                    br=get_birthrate(m.state) 
                    ci = clone(pop,m.state,br)
                    s.popDict[m.state]=ci
                else
                    push!(s.popDict[m.state].pop,c)
                end
            # else... ASSUME missegregations from stochastic compartment into det. compartment
            # are negligible (!?)              
            end
        end
        if length(s.popDict[k].pop) < 1
            delete!(s.popDict,k)            
        end
    end
    return(s)
end
# perform cell death for entire stochastic model
function die(s::stochasticCompartment)
    for k in keys(s.popDict)
        s.popDict[k] = die(s.popDict[k],s.dt,s.deathRate,s.uniform01)
        if length(s.popDict[k].pop) < 1
            delete!(s.popDict,k)            
        end
    end
    return(s)
end

# return grid with locations of all cells in stochastic compartment
function spatial_summary(s::stochasticCompartment)
    u=zeros(s.Np...)
    for cl in values(s.popDict)
        for c in cl.pop
            u[ceil(Int64,c.x),ceil(Int64,c.x)] = u[ceil(Int64,c.x),ceil(Int64,c.x)]+1
        end
    end
    return u
end

# introduce new cells into stochastic model from deterministic compartment
function getDaughters(stoch::stochasticCompartment, E0, s0, sIndex,birthRates)
    for i in 1:(size(s0)[1])
        br = birthRates[i]
        for j in 1:(size(s0)[2])
            for k in 1:(size(s0)[3])
                pBirth = br*stoch.dt*s0[i,j,k]*E0[j,k]/(stoch.ϕ + E0[j,k]) # expected number of births per site
                pMis = pBirth*stoch.misrate # expected number of missegregations per site
                p = Poisson(pMis)
                Np=rand(p,1)[1]
                if Np>0
                    for n in 1:Np 
                        k0=karyotype(deepcopy(sIndex[i]),true)
                        k1,k2=missegregate(k0, stoch.maxChrom, stoch.chromosome_selector)
                        for ki in [k1,k2]        
                            if ki.viable
                                # create a new cell
                                c=cell((j-rand()),(k-rand()))
                                # check if cell karyotype exists in clones
                                if (ki.state in keys(stoch.popDict)) == false
                                #popDict[k.state]=1
                                    pop = MutableLinkedList{cell}()
                                    push!(pop, c)
                                    br=get_birthrate(ki.state) # will look this up later
                                    ci = clone(pop,ki.state,br)
                                    stoch.popDict[ki.state]=ci
                                else
                                    push!(stoch.popDict[ki.state].pop,c)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return(stoch)
end

# get x,y & cn of each cell and return in a dataframe
function summarise_clone(ci::clone)
    df = repeat(ci.cn',outer=(length(ci.pop),1))
    xy = hcat(map(c->[c.x,c.y],ci.pop)...)'
    return DataFrame(hcat(df,xy),:auto)
end

## should take as arguments the energy array, the cell array, and stochasticCOmpartmentobject
## then return updated versions of each
function stochastic_step(s::stochasticCompartment, E0, s0, info)
    @unpack (sIndex,birthRates) = info
    dE0 = grad(E0,s.Np)
    #unclear what is best order to do these steps
    s=die(s)
    s=migrate(s,dE0,E0)
    s=replicate(s,E0,dE0,sIndex)
    s=getDaughters(s, E0, s0, sIndex,birthRates)
    u=spatial_summary(s)
    return s, u
end

## wrapper function to setup the stochastic model
## pass to the functions variables that remain constant 
## throughout simulation etc. Should return a stochasticCompartment object
function setup_stochastic(params)

    @unpack (Np,ϕ, nChrom,maxChrom,dt,misrate,deathRate,Γ,χ,Ξ) = params
    chromosome_selector = DiscreteUniform(1,nChrom) # random integer 1:nchrom
    uniform01 = Uniform(0,1) # random continous uniform dist.
    dnorm = Normal(0,dt*Γ) #normal distribution for cell migration

    popDict = Dict{Array,clone}()
    s = stochasticCompartment(popDict,Np,ϕ,χ,Ξ, maxChrom,dt,misrate,deathRate,
    chromosome_selector,uniform01,dnorm)
    return(s)
end

# implementation of wrapped boundary coords
function wrap(index,Npi)
    if index>Npi
        return index-Npi
    end
    if index<1
        return index+Npi
    end
    return index
end

#get the gradient of a 2d array
function grad(R2,Np)
    dR2 = zeros(2,size(R2)...)
    for i in 1:size(R2)[1]
        for j in 1:size(R2)[2]
            #1st order central diff
            dR2[1,i,j]=(R2[wrap(i+1,Np[1]),j]-R2[wrap(i-1,Np[1]),j])/(2*Np[1])
            dR2[2,i,j]=(R2[i,wrap(j+1,Np[2])]-R2[i,wrap(j-1,Np[2])])/(2*Np[2])
        end
    end
    return dR2
end


