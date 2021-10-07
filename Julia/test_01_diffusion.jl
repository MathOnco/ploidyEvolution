using DifferentialEquations,ComponentArrays

## test cell diffusion in PDE model
include("stochastic.jl")
include("ploidyMovement.jl")

Γtest=500
tspan=(0.0,10.0)
Np=[20,20]
Lp=[1000,1000]

boundary_Dict = Dict{Tuple{Int,Int},Any}()
domain_Dict = Dict{Tuple{Int,Int},Int}()

dp = Lp./(Np .- 1)
starting_copy_number=[1,1,1,1,1]
sIndex=[starting_copy_number]
birthRates = [0.37]

odePars = (nChrom=5,
misRate=0.0,
deathRate=0.0,
Γ=Γtest,Γₑ=0.0,ϕ=0.1,Ξ=0.1,χ=0.0,δ=0.0,Np=Np,dp=dp,k=0.0,E_vessel=0.0,
domain_Dict=domain_Dict,
boundary_Dict=boundary_Dict,
max_cell_cycle_duration=5,
max_birthRate=0.37,
debugging=0)
M=zeros(1,1)
s0=zeros(1,Np...)
s0[1,ceil(Int64,Np[1]/2),ceil(Int64,Np[2]/2)]=100000
u0=s0[1,:,:]
open("test_output/01_diffusion/u0.csv","w") do io
    writedlm(io,u0,',')
end

E0=zeros(Np...)
u_s=zeros(Np...)
u0 = ComponentArray(s=s0,E=E0)
odePars=(M=M,u_s=u_s,birthRates=birthRates,sIndex=sIndex,odePars...)

prob = ODEProblem(ploidy_hybrid_model,u0,tspan,odePars)
sol = solve(prob,VCABM(),abstol=1e-8,reltol=1e-5)
#println("Stochastic step...")
s=sol[length(sol)].s

open("test_output/01_diffusion/pde.csv","w") do io
    for i in 1:length(sIndex)
        z = s[1,:,:]
        writedlm(io,z,',')
    end
end

## test cell diffusion in stochastic model

pop = MutableLinkedList{cell}()
br=0.0
dt=0.1
χ=0.0
Ξ=0.1
Γ=Γtest/max(dp...)^2
dnorm = Normal(0,sqrt(2*Γ*dt)) 
for i in 1:100000
    c = cell(ceil(Int64,Np[1]/2)-rand(),ceil(Int64,Np[2]/2)-rand())
    push!(pop, c)
end
ci = clone(pop,starting_copy_number,br)

u=zeros(Np...)
dE0=zeros(2,Np...)

for i in 1:100
    migrate(ci,Np,dnorm,dE0,E0,χ,Ξ, domain_Dict)
end

for cl in ci.pop
    u[ceil(Int64,cl.x),ceil(Int64,cl.y)]=u[ceil(Int64,cl.x),ceil(Int64,cl.y)]+1
end

open("test_output/01_diffusion/stochastic.csv","w") do io
    writedlm(io,u,',')
end