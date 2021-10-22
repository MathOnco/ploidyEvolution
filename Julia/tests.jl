function load_vessels(Lp,Np,vfile)
    df = CSV.read(vfile,DataFrame)
    x = range(0,Lp[1],length=Np[1])
    y = range(0,Lp[2],length=Np[2])
    dp = Lp./(Np .- 1)

    domain_Dict = Dict{Tuple{Int,Int},Int}() # 1 = boundary blood vessel, 2 = inside blood vessel

    for i = 1 : Np[1]
        xdist = (x[i] .- df[!,"Centroid X µm"]).^2
        for j = 1 : Np[2]
            ydist = (y[j] .- df[!,"Centroid Y µm"]).^2
            if minimum(xdist .+ ydist) < sum(dp.^2)
                domain_Dict[(i,j)]=1
            end
        end
    end

    # But first we find whether the point is interior
    for k in keys(domain_Dict)
        # Vessel interior point
        if issubset(grid_neighbours(k,Np[1]),keys(domain_Dict))
            domain_Dict[k] = 2
        end
    end

    for (k,v) in domain_Dict
        if v == 1
            delete!(domain_Dict,k)
        end
    end

    vesselkeys=keys(domain_Dict)

    # Fill dict with info about the normal vector
    
    interface_i=[]
    interface_o=[]
    vessels=[]
    for v in vesselkeys
        vessels=[vessels...,v]
        for n in grid_neighbours(v,Np[1])
            if n in vesselkeys
            else
                interface_o=[interface_o...,n]
                interface_i=[interface_i...,v]
            end
        end

    end

    interfaces=hcat(interface_o,interface_i)

    return vessels,interfaces

end

function taxis_ij(coord,Npi,E,s,s0,Ξ,dpi,χ)
    #chemotaxis across boundary between (i,j) and (i+1,j)
    dE = chemotaxis_form(E[wrap(coord[1]+1,Npi),coord[2]],Ξ)-chemotaxis_form(E[coord],Ξ)
    if dE>0
        qi= χ*dE*map(si->Skplus(s[si,wrap(coord[1]-1,Npi),coord[2]],s[si,coord],s[si,wrap(coord[1]+1,Npi),coord[2]]),1:size(s,1))/dpi^2
        s0[:,coord].-=qi
        s0[:,wrap(coord[1]+1,Npi),coord[2]].+=qi
    else
        qi=χ*dE*map(si->Skminus(s[si,coord],s[si,wrap(coord[1]+1,Npi),coord[2]],s[si,wrap(coord[1]+2,Npi),coord[2]]),1:size(s,1))/dpi^2
        s0[:,coord]-=qi
        s0[:,wrap(coord[1]+1,Npi),coord[2]]+=qi
    end
    #chemotaxis across boundary between (i,j) and (i,j+1)
    dE = chemotaxis_form(E[coord[1],wrap(coord[2]+1,Npi)],Ξ)-chemotaxis_form(E[coord],Ξ)
    if dE>0
        qi= χ*dE*map(si->Skplus(s[si,coord[1],wrap(coord[2]-1,Npi)],s[si,coord],s[si,coord[1],wrap(coord[2]+1,Npi)]),1:size(s,1))/dpi^2
        s0[:,coord].-=qi
        s0[:,coord[1],wrap(coord[2]+1,Npi)].+=qi
    else
        qi=χ*dE*map(si->Skminus(s[si,coord],s[si,coord[1],wrap(coord[2]+1,Npi)],s[si,coord[1],wrap(coord[2]+2,Npi)]),1:size(s,1))/dpi^2
        s0[:,coord]-=qi
        s0[:,coord[1],wrap(coord[2]+1,Npi)]+=qi
    end

end

birth_ij = function(coord,s,s0,E,M,birthrates,ϕ,max_cell_cycle_duration)
    br = energy_constrained_birth(E[coord],ϕ)*birthrates
    br = br.*(br.>log(2)/max_cell_cycle_duration)
    s0[:,coord].+=M*(s[:,coord].*br)
end

diffu_ij = function(coord,s,s0,Γ,dpi,Npi)
    s0[:,coord].+=hcat(map(x->Γ*s[:,x...],grid_neighbours(coord,Npi))...)*ones(4)/dpi^2-
    Γ*4*s[:,coord]/dpi^2
end

diffE_ij = function(coord,E,dE,Γₑ,dpi,Npi)
    dE[coord].+=sum(map(x->Γₑ*E[x...],grid_neighbours(coord,Npi)))/dpi^2-
    Γₑ*4*E[coord]/dpi^2
end

function noflux(exterior,interior,E,dE,s,ds,dpi,Npi,Γ,χ,Ξ,Γₑ,k,E_vessel)

    ##vessel behaviour
   # dE[exterior...]+=Γₑ*E[exterior...]/dpi^2-Γₑ*E_vessel/dpi^2 #+
    #Γₑ*k*(E_vessel-E[exterior...])/dpi^2

    #noflux in diffusion
    ## we can simply "give back" 25% of the efflux per bounding coord:
    ds[:,exterior...]=ds[:,exterior...]+Γ*s[:,exterior...]/dpi^2 

    ##noflux in chemotaxis
    dE = chemotaxis_form(E[interior...],Ξ)-chemotaxis_form(E[exterior...],Ξ)
    if dE<0
        println("Assumption that vessel has more energy than tissue is violated!")
    else
        np1 = (wrap(2*exterior[1]-interior[1],Npi),wrap(2*exterior[2]-interior[2],Npi))

        if sum(exterior)>sum(interior) ## we are correcting migration across the vessel in a negative direction
            ds[:,exterior...]+=χ*dE*map(si->Skminus(s[si,interior...],s[si,exterior...],s[si,np1...]),1:size(s,1))/dpi^2
        else ## we are correcting migration across the vessel in a positive direction
            ds[:,exterior...]+= χ*dE*map(si->Skplus(s[si,np1...],s[si,exterior...],s[si,interior...]),1:size(s,1))/dpi^2
        end

    end

end

function interiors(coord,dE,ds)
    dE[coord...]=0
    ds[:,coord...].=0
end

function consumption(coord,dE,E,s,u_s,δ)
    dE[coord] -= δ*E[coord]*(sum(s[:,coord])+u_s[coord])
end


function ploidy_hybrid_model(du,u,pars,t)
    #println(t)
	# Grab the parameters
	@unpack (domain_Dict, boundary_Dict,M,u_s,birthRates,sIndex,nChrom,misRate,deathRate,Γ,Γₑ,ϕ,Ξ,χ,δ,Np,dp,k,E_vessel,
	vessels,interfaces,max_cell_cycle_duration,max_birthRate,debugging) = pars

    nChrom = length(sIndex[1])

    s, E = u.s, u.E
    Npi = Np[1]
    dpi=dp[1]

    coordsList = [1:num_points for num_points in Np]
	# dp = (Np .- 1).^(-1)
	sum_dp_invsq = sum(dp.^(-2))

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

        diffusion = Γₑ*(sum((E_plus[i] + E_minus[i])/dp[i]^2 for i in 1:length(dp)) - 2.0*sum_dp_invsq*E_focal)

		# Consumption by the cells at the grid point designated by coord.
		consumption = δ*E_focal*(sum(s[:,coord...])+u_s[coord...])

		# update the energy compartment
		du.E[coord...] = diffusion - consumption

    end


    
	map(i->taxis_ij(i,Npi,E,s,du.s,Ξ,dpi,χ),CartesianIndices(E))
    map(x->diffu_ij(x,s,du.s,Γ,dpi,Npi),CartesianIndices(E))
    map(x->birth_ij(x,s,du.s,E,M,birthRates,ϕ,max_cell_cycle_duration),CartesianIndices(E))
    du.s -= deathRate*s #death
    #map(x->consumption(x,du.E,E,s,u_s,δ),CartesianIndices(E))
    #map(x->diffE_ij(x,E,du.E,Γₑ,dpi,Npi),CartesianIndices(E))
    map(i->noflux(interfaces[i,1],interfaces[i,2],E,du.E,s,du.s,dpi,Npi,Γ,χ,Ξ,Γₑ,k,E_vessel),1:size(interfaces,1))
    map(x->interiors(x,du.E,du.s),vessels)


end




