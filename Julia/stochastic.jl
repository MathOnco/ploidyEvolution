using DataStructures, Distributions, CSV, DataFrames, Random, Parameters

## TO DO
# make move() failsafe
# input a threshold max number of cells
# increase compatibility with ODE model (e.g. polyharmonicSpline birth rates; use chemotaxis_form() fn's; & etc)
# modify birth_rate() to allow testing with different maxChrom 
# need to pass the boundary sites so that cells cannot migrate into vessels
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
function get_birthrate(cn,interp,option=1)
    if option==1
        return max(PolyharmonicInterpolation.polyharmonicSpline(interp,cn')[1],0.0)
    end
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
function move(c::cell,Np,dnorm, dE0,E0,χ,Ξ,domain_Dict)
    if get(domain_Dict,(ceil(Int64,c.x),ceil(Int64,c.y)),0) == 2
        println("$c in vessel! move()")
    end
    
    ix = ceil(Int64,c.x)
    iy = ceil(Int64,c.y)
    
    new_x = 0.
    new_y = 0.
    # do not allow cells to move into vessels
    # question is how to achieve this... i.e. what is desired behaviour? A "fair" thing to do
    # would be to keep generating random sets of new coordinates until a set is drawn that is not inside a vessel.
    # however, in practise with chemotaxis in play this could create situations where cells near vessels virtually always want to go inside a vessel
    # and may even end up frequently skipping "over" vessels. Therefore strategy will jsut be to have one attempt at coords
    # and if they end up inside a vessel, cell just retains its old coords.  
    
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
    
    if get(domain_Dict,(ceil(Int64,new_x),ceil(Int64,new_y)),0) != 2
        c.x = new_x
        c.y = new_y
    end

    
    return c
end

#store all the cells for a given clone
mutable struct clone
    pop
    cn
    birthrate
end

# wrapper function for move()
function migrate(cl::clone,Np,dnorm,dE0,E0,χ,Ξ, domain_Dict)
    cl.pop = map(x->move(x,Np,dnorm,dE0,E0,χ,Ξ,domain_Dict),cl.pop)
    return(cl)
end

# loop thru all cells per clone, and perform divisions
function replicate(cl::clone, E0, ϕ, dt, uniform01,misrate,maxChrom,chromosome_selector,Np,dnorm,dE0,χ,Ξ,max_cell_cycle_duration,domain_Dict)
    # for missegregations we will take a strategy of returning all the new viable karyotypes
    # and then allowing the stochasticCompartment object to deal with them
    missegs = MutableLinkedList{cell_allinfo}()
    for i in length(cl.pop):-1:1
        energy = E0[ceil(Int64,cl.pop[i].x),ceil(Int64,cl.pop[i].y)]
        #need to replace with call to energy function
        energy_modified_birth_rate =cl.birthrate*energy/(ϕ + energy)
        if log(2)/max_cell_cycle_duration < energy_modified_birth_rate  && rand(uniform01,1)[1] < energy_modified_birth_rate*dt # division occurs
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
                new_c=move(new_c,Np,dnorm,dE0,E0,χ,Ξ,domain_Dict)
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
    ##
    interp
    max_cell_cycle_duration
    domain_Dict
    stochasticThreshold
end

# perform migration for entire stochastic model
function migrate(s::stochasticCompartment,dE0,E0)
    for k in keys(s.popDict)
        s.popDict[k] = migrate(s.popDict[k],s.Np,s.dnorm,dE0,E0,s.χ,s.Ξ,s.domain_Dict)
    end
    return(s)
end
# perform cell division for entire stochastic model
function replicate(s::stochasticCompartment, E0,dE0,sIndex,s0)
    for k in keys(s.popDict)
        s.popDict[k],missegs = replicate(s.popDict[k],E0, s.ϕ, s.dt, s.uniform01,
        s.misrate,s.maxChrom,s.chromosome_selector,s.Np,s.dnorm,dE0,s.χ,s.Ξ,s.max_cell_cycle_duration,s.domain_Dict)
        for m in missegs
            c = cell(m.x,m.y)
            c=move(c,s.Np,s.dnorm,dE0,E0,s.χ,s.Ξ,s.domain_Dict)
            if (m.state in sIndex)==false # check if already in det. comp.
                if (m.state in keys(s.popDict)) == false 
                    pop = MutableLinkedList{cell}()
                    push!(pop, c)
                    br=get_birthrate(m.state,s.interp) 
                    ci = clone(pop,m.state,br)
                    s.popDict[m.state]=ci
                else
                    push!(s.popDict[m.state].pop,c)
                end
            else
                si=findall(==(m.state),sIndex)
                s0[si,ceil(Int64,c.x),ceil(Int64,c.y)].+=1
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
            u[ceil(Int64,c.x),ceil(Int64,c.y)] = u[ceil(Int64,c.x),ceil(Int64,c.y)]+1
        end
    end
    return u
end


function graduate_clones(s::stochasticCompartment)
    #currently just allowing one clone to graduate per timestep
    for k in keys(s.popDict)
        if length(s.popDict[k].pop) > s.stochasticThreshold 
            clones = k 
            birthRate = s.popDict[k].birthrate
            pos=zeros(1,s.Np...)
            for c in s.popDict[k].pop
                pos[1,ceil(Int64,c.x),ceil(Int64,c.y)]=pos[1,ceil(Int64,c.x),ceil(Int64,c.y)]+1
            end
            delete!(s.popDict,k)
            return clones, birthRate, pos
        end
    end
    clones=0
    birthRate=0
    pos=0
    return clones,birthRate, pos
end

# place a clone population defined on a grid into the stochastic model
function addClone(s::stochasticCompartment, ui, cn, br)
    if cn in keys(s.popDict) 
        error("population already in stochastic compartment!")
    end
    pop = MutableLinkedList{cell}()
    ci = clone(pop,cn,br)
    for i in (1:size(ui)[1])
        for j in (1:size(ui)[2])
            # different possibilities exist for how to transfer deterministic->stochastic
            # here we choose to round the deterministic population to the nearest Int per grid site and emplace into stoch. model
            # but another alternative could be to treat det. pop. as a probability distribution and randomly realise from that dist
            Ncells = round(ui[i,j])
            if Ncells>0
                for k in 1:Ncells
                    c=cell((i-rand()),(j-rand()))
                    push!(ci.pop, c)
                end
            end
        end
    end
    if  length(ci.pop)<1 
        error("tried to add population of zero cells to stochastic compartment!")
    end
    s.popDict[cn]=ci
    return(s)
end

# introduce new cells into stochastic model from deterministic compartment
function getDaughters(stoch::stochasticCompartment, E0, s0, sIndex,birthRates)
    for i in 1:(size(s0)[1])
        br0 = birthRates[i]
        for j in 1:(size(s0)[2])
            for k in 1:(size(s0)[3])
                energy_modified_birth_rate = br0*E0[j,k]/(stoch.ϕ + E0[j,k])
                if log(2)/stoch.max_cell_cycle_duration < energy_modified_birth_rate
                    pBirth = energy_modified_birth_rate*stoch.dt*s0[i,j,k] # expected number of births per site
                    pMis = pBirth*stoch.misrate # expected number of missegregations per site
                    p = Poisson(pMis)
                    Np=rand(p,1)[1]
                    if Np>0
                        for n in 1:Np 
                            k0=karyotype(deepcopy(sIndex[i]),true)
                            k1,k2=missegregate(k0, stoch.maxChrom, stoch.chromosome_selector)
                            for ki in [k1,k2]        
                                if ki.viable && (ki.state in sIndex) == false
                                    # create a new cell
                                    c=cell((j-rand()),(k-rand()))
                                    if get(stoch.domain_Dict,(ceil(Int64,c.x),ceil(Int64,c.y)),0) == 2
                                        println("$c in vessel! daughters()... (j,k)")
                                    end
                                    # check if cell karyotype exists in clones
                                    if (ki.state in keys(stoch.popDict)) == false 
                                    #popDict[k.state]=1
                                        pop = MutableLinkedList{cell}()
                                        push!(pop, c)
                                        br=get_birthrate(ki.state,stoch.interp) # will look this up later
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
    end
    return(stoch)
end

# get x,y & cn of each cell and return in a dataframe
function summarise_clone(ci::clone)
    df = repeat(ci.cn',outer=(length(ci.pop),1))
    xy = hcat(map(c->[c.x,c.y],ci.pop)...)'
    return DataFrame(hcat(df,xy),:auto)
end

# get x,y & cn of each cell and return in a dataframe
function count_clones(s::stochasticCompartment)
    df = map(ci->hcat(ci.cn...,ci.birthrate,length(ci.pop)),values(s.popDict))
    df = reduce(vcat,df)
    return df
end

# overload function def for use also by pde model
function count_clones(s,sIndex,birthRates)
    df = map(i->hcat(sIndex[i]...,birthRates[i],sum(s[i,:,:])),1:length(sIndex))
    df = reduce(vcat,df)
    return df
end

## should take as arguments the energy array, the cell array, and stochasticCOmpartmentobject
## then return updated versions of each
function stochastic_step(s::stochasticCompartment, E0, s0, info)
    @unpack (sIndex,birthRates) = info
    dE0 = grad(E0,s.Np)
    #unclear what is best order to do these steps
    s=die(s)
    s=migrate(s,dE0,E0)
    s=replicate(s,E0,dE0,sIndex,s0)
    s=getDaughters(s, E0, s0, sIndex,birthRates)
    u=spatial_summary(s)
    new_clone, new_birthRate, new_pos = graduate_clones(s)
    return s, u, new_clone, new_birthRate, new_pos
end

## wrapper function to setup the stochastic model
## pass to the functions variables that remain constant 
## throughout simulation etc. Should return a stochasticCompartment object
function setup_stochastic(params)

    @unpack (Lp,Np,ϕ, nChrom,maxChrom,dt,misrate,deathRate,Γ,χ,Ξ,interp,max_cell_cycle_duration, domain_Dict, stochasticThreshold) = params
    
    # needs fixed. We need to trace back these parms and determine their units
    dp=Lp./(Np .- 1)
    Γ=Γ/max(dp...)^2
    χ=χ/max(dp...)^2
    
    chromosome_selector = DiscreteUniform(1,nChrom) # random integer 1:nchrom
    uniform01 = Uniform(0,1) # random continous uniform dist.
    dnorm = Normal(0,dt*Γ) #normal distribution for cell migration

    popDict = Dict{Array,clone}()
    s = stochasticCompartment(popDict,Np,ϕ,χ,Ξ, maxChrom,dt,misrate,deathRate,
    chromosome_selector,uniform01,dnorm,interp,max_cell_cycle_duration, domain_Dict, stochasticThreshold)
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

#implement wrapped boundary conditions on grid
function grid_neighbours(coord,Npi)
    nx1= (wrap(coord[1]-1,Npi),coord[2])
    nx2= (wrap(coord[1]+1,Npi),coord[2])
    ny1= (coord[1],wrap(coord[2]-1,Npi))
    ny2= (coord[1],wrap(coord[2]+1,Npi))

    neighbours = (nx1,nx2,ny1,ny2)
    return neighbours
end

function ploidy_hybrid_model(du,u,pars,t)

	# Grab the parameters
	@unpack (M,u_s,birthRates,sIndex,nChrom,chromArray,misRate,deathRate,Γ,Γₑ,ϕ,Ξ,χ,δ,Np,dp,k,E_vessel,
	domain_Dict,boundary_Dict,max_cell_cycle_duration,max_birthRate,debugging) = pars

    nChrom = length(sIndex[1])

	# Convert chrom array to integer
	intChromArray = [Int(first(x)):Int(step(x)):Int(last(x)) for x in chromArray]
	s, E = u.s, u.E

	coordsList = [1:num_points for num_points in Np]
	# dp = (Np .- 1).^(-1)
	sum_dp_invsq = sum(dp.^(-2))

	# cflow = 0.0
	# dflow = 0.0

	# Iterate over each coordinate
	for coord in Iterators.product((1:length(coords) for coords in coordsList)...)

		# Get whether we are interior (0), at boundary of vessel (1) or inside vessel (2)
		domain_info = get(domain_Dict,coord,0)

		if domain_info == 2 # Inside blood vessel we make no changes
			du.s[:,coord...] .= 0.0
			du.E[coord...] = 0.0
			continue
		end

		# Gets indicies used to calculate spatial impact (assume periodic)
		cardinal_point_info = get_cardinal_gridpoints([elem for elem in coord],Np)

		# Gets the points to be used for spatial grid resolution
		E_cardinal = [E[info...] for info in cardinal_point_info]
		E_minus,E_plus = E_cardinal[end-1:-2:1],E_cardinal[end:-2:1]
		E_focal = E[coord...]

		if domain_info == 1	# boundary of vessel
			normal_vector = boundary_Dict[coord] 	# Grab normal vector
			for i = 1 : length(dp)
				if normal_vector[i] < 0 # The rightward point is in the vessel so we replace it with the B.C.
					E_plus[i] = E_minus[i] - 2.0*dp[i]*k*(E_focal - E_vessel)
				elseif normal_vector[i] > 0 # The leftward point is in the vessel
					E_minus[i] =  E_plus[i] - 2.0*dp[i]*k*(E_focal - E_vessel)
				end
			end
		end

		# Iterate over each compartment
		for si in 1:length(sIndex)

			# Initialize inflow to 0
			inflow = 0.0

			# copy state to local variables for brevity
			s_focal = s[si...,coord...]
			s_cardinal = [s[si...,info...] for info in cardinal_point_info]
			s_minus,s_plus = s_cardinal[end-1:-2:1],s_cardinal[end:-2:1]

			# Checking the boundary
			if domain_info == 1	# boundary of vessel
				normal_vector = boundary_Dict[coord] 	# Grab normal vector
				for i = 1 : length(dp)
					if normal_vector[i] < 0 # The rightward point is in the vessel so we replace it with the B.C.
						s_plus[i] = s_minus[i] + χ/Γ*s_focal*(chemotaxis_form(E_plus[i],Ξ) - 
						chemotaxis_form(E_minus[i],Ξ))
					elseif normal_vector[i] > 0 # The leftward point is in the vessel
						s_minus[i] = s_plus[i] - χ/Γ*s_focal*(chemotaxis_form(E_plus[i],Ξ) - 
						chemotaxis_form(E_minus[i],Ξ))
					end
				end
			end

			energy_constraint = energy_constrained_birth(E_focal,ϕ)
			
            for sj in 1:length(sIndex)

                s_parent = s[sj,coord...]

                # Get max birth rate
                birthRate = M[sj,si]

                if s_parent == 0.0 || birthRate == 0.0
                    continue
                end

                # If birth rate is below a threshold determined by amount spent in G1 (default is 100 hours)
                # then we assume that birth does not occur...
                energy_modified_birth_rate = birthRate*energy_constraint
                if energy_modified_birth_rate > log(2)/max_cell_cycle_duration
                    # Get flow rate from parentCN -> focalCN
                    flowRate = misRate/nChrom
                    inflow += energy_modified_birth_rate*flowRate*s_parent
                else
                    # println("At energy level $(E_focal), birth rate of $(parentCN) is too low...")
                end

            end


            if s_focal > 0
                # Grab the focal cells (max?) birth rate
                birthRate = birthRates[si]
                energy_modified_birth_rate = birthRate*energy_constraint

                if energy_modified_birth_rate > log(2)/max_cell_cycle_duration
                    # Add it to the inflow to the focal compartment
                    inflow += energy_modified_birth_rate*(1.0 - 2.0*misRate)*s_focal
                else
                    # println("At energy level $(E_focal), birth rate of $(focalCN) is too low...")
                end
            end

			# Death of the focal compartment
			outflow = deathRate*s_focal

			# u_xx -> (u[i+1] + u[i-1] - 2u[i])/dx^2 
			# random migration
			diffusion = Γ*(sum((s_plus[i] + s_minus[i])/dp[i]^2 for i in 1:length(dp)) 
			- 2.0*sum_dp_invsq*s_focal)

			# chemotaxis
			chemotaxis = χ*(
				s_focal*(
				sum(
					(chemotaxis_form(E_plus[i],Ξ) + chemotaxis_form(E_minus[i],Ξ))/dp[i]^2
					for i in 1 : length(dp)
					)
				- 2.0*sum_dp_invsq*chemotaxis_form(E_focal,Ξ)) + 
				sum(
					(s_plus[i] - s_minus[i])*(chemotaxis_form(E_plus[i],Ξ) - chemotaxis_form(E_minus[i],Ξ))/dp[i]^2
					for i in 1 : length(dp)))/4.0

			# Update the RHS dependent on where we are in the domain
			du.s[si,coord...] = inflow - outflow + diffusion - chemotaxis

			# dflow += diffusion
			# cflow += chemotaxis

		end

		# Diffusion of energy molecule through the tissue
		diffusion = Γₑ*(sum((E_plus[i] + E_minus[i])/dp[i]^2 for i in 1:length(dp)) - 2.0*sum_dp_invsq*E_focal)

		# Consumption by the cells at the grid point designated by coord.
		consumption = δ*E_focal*(sum(s[:,coord...])+u_s[coord...])

		# update the energy compartment
		du.E[coord...] = diffusion - consumption

		# dflow += diffusion
		
	end

	# @show dflow,cflow

end

#calculate transition matrix for missegregations between states in pde model
function transition_matrix(sIndex,birthRates)

    M = zeros(length(sIndex),length(sIndex))
    for i in 1:length(sIndex)
        for j in 1:length(sIndex)
            if sum(map(abs,(sIndex[i]-sIndex[j])))==1
                M[i,j] = birthRates[i]
            end
        end
    end
    return M
end

function run_hybrid_step(i,odePars,u,tspan,s_s,saveat,cbset,sIndex,birthRates,M,u_s)
	#M as well will need to be pulled out of odePars
	#println("PDE step...")
	odePars=(M=M,u_s=u_s,birthRates=birthRates,sIndex=sIndex,odePars...)
	prob = ODEProblem(ploidy_hybrid_model,u,tspan,odePars)
	sol = solve(prob,VCABM(),abstol=1e-8,reltol=1e-5,saveat=saveat,callback=cbset)
	#println("Stochastic step...")
	u.E=sol[length(sol)].E
	u.s=sol[length(sol)].s


	# check if any PDE compartments have become too small
	compartment_sizes = map(x->sum(u.s[x,:,:]),1:length(sIndex))
	relegation_index=findmin(compartment_sizes)

	stochInfo=(sIndex=sIndex,birthRates=birthRates)
	s_s, u_s, new_clone, newBirthRate, new_pos =stochastic_step(s_s,sol[length(sol)].E,sol[length(sol)].s,stochInfo)
	s_new=u.s
	append_state=false
	if (new_clone==0)==false
		append_state=true
		sIndex = (sIndex..., new_clone)
		birthRates=(birthRates...,newBirthRate)
		s_new=cat(u.s,new_pos,dims=1)
	end
	return s_new, sIndex, birthRates, relegation_index, u_s, append_state
end