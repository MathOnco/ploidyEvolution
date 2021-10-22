using DifferentialEquations,ComponentArrays, Plots

## test cell diffusion in PDE model
include("stochastic.jl")
include("ploidyMovement.jl")

Γtest=100
χtest=1000
tspan=(0.0,100.0)
Np=[30,30]
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
M=zeros(1,1)
s0=zeros(1,Np...)
s0[:,ceil(Int64,Np[1]/2),ceil(Int64,Np[1]/4)].=100000
E0=zeros(Np...)

for i in 1:Np[1]
    pp = Lp[1]*(i)/Np[1]/100
    E0[:,i].= exp(pp)-0.02
end

u_s=zeros(Np...)
u0 = ComponentArray(s=s0,E=E0)
odePars=(M=M,u_s=u_s,birthRates=birthRates,sIndex=sIndex,odePars...)

prob = ODEProblem(ploidy_hybrid_model,u0,tspan,odePars)
sol = solve(prob,AutoTsit5(Rosenbrock23()))
#println("Stochastic step...")
for t in 0:10:round(Int64,tspan[2])
    s=sol(t).s
    fn="test_output/03_taxis/linear/"*string(t,pad=3)*"cells.png"
    xt = heatmap(sol(t).s[1,:,:])
    display(xt)
    png(fn) 
    open("test_output/03_taxis/linear/"*string(t,pad=3)*"pde.csv","w") do io
        for i in 1:length(sIndex)
            z = s[1,:,:]
            writedlm(io,z,',')
        end
    end
    

end

## test cell diffusion in stochastic model

pop = MutableLinkedList{cell}()
br=0.0
dt=0.1
χ=dt*χtest/max(dp...)

Γ=Γtest/max(dp...)^2
dnorm = Normal(0,sqrt(2*Γ*dt)) 
for i in 1:100000
    c = cell(ceil(Int64,Np[1]/2)-rand(),ceil(Int64,Np[1]/4)-rand())
    
    push!(pop, c)
end
ci = clone(pop,starting_copy_number,br)

u=zeros(Np...)


dE0=map(x->chemotaxis_form(x,Ξ),E0)
dE0 = grad(dE0,Np,dp)

for i in 1:ceil(Int64,tspan[2]/dt)
    migrate(ci,Np,dnorm,dE0,E0,χ,Ξ, domain_Dict)
    if i%100==0
        u=zeros(Np...)
        for cl in ci.pop
            u[ceil(Int64,cl.x),ceil(Int64,cl.y)]=u[ceil(Int64,cl.x),ceil(Int64,cl.y)]+1
        end
        savestamp=round(Int64,i/10)
        open("test_output/03_taxis/linear/"*string(savestamp,pad=3)*"stoch.csv","w") do io
            writedlm(io,u,',')
        end
    end
end


open("test_output/03_taxis/linear/E0.csv","w") do io
    writedlm(io,E0,',')
end