using DifferentialEquations,ComponentArrays, Plots, Profile, BenchmarkTools
# purpose of this file is to test  the execution speed of the main 
# ploidy PDE function. 

include("ploidyMovement.jl")


Γtest=100
χtest=1000
tspan=(0.0,50.0)
Np=[50,50]
Lp=[1000,1000]
Ξ=0.02

boundary_Dict = Dict{Tuple{Int,Int},Any}()
domain_Dict = Dict{Tuple{Int,Int},Int}()

dp = Lp./(Np .- 1)
pp=0:(Lp[1]/(Np[1]-1)):Lp[1]
starting_copy_number=[1,1,1,1,1]
sIndex=[starting_copy_number]
birthRates = [0]

odePars = (nChrom=5,
misRate=0.0,
deathRate=0.0,
Γ=Γtest,Γₑ=0.0,ϕ=0.02,Ξ=Ξ,χ=χtest,δ=0.0,Np=Np,dp=dp,k=0.0,E_vessel=0.0,
domain_Dict=domain_Dict,
boundary_Dict=boundary_Dict,
max_cell_cycle_duration=5,
max_birthRate=0,
debugging=0)
M=zeros(2,2)
s0=rand(2,Np...)
E0=rand(Np...)
ds=rand(2,Np...)
dE=rand(Np...)


u_s=rand(Np...)
u0 = ComponentArray(s=s0,E=E0)
du=ComponentArray(s=ds,E=dE)

odePars=(M=M,u_s=u_s,birthRates=birthRates,sIndex=sIndex,odePars...)

function run10x(du,u0,odePars)
    for i in 1:10
        #println(i)
        ploidy_hybrid_model(du,u0,odePars,0)
    end
end


@btime run10x(du,u0,odePars)
function run10x_o(du,u0,odePars)
    for i in 1:10
        ploidy_hybrid_model_old(du,u0,odePars,0)
    end
end


@btime run10x_o(du,u0,odePars)

