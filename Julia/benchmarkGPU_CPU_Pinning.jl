using Distributed
using BenchmarkTools
using CUDA
using StatsBase
using Random


function test(ttype::String, n::Int64)
    B_ = ones(Float32,n);
    B = cu(B_);
    A = rand(Float32, length(B_));

    for i in 1:20
        if ttype == "cpu"
			o = B_.^3 + B_ - B_
        elseif ttype == "gpu"
            o = B.^3  + B - B
        elseif ttype == "gpu broadcast"
            B[1:100] .= 2
            o = B.^3  + B - B
        elseif ttype == "gpu copy"
            B = cu(A)
            o = B.^3  + B - B
        elseif ttype == "gpu pin"
            A = Mem.pin(Array{eltype(A)}(undef, size(A)));
            Ap = pointer(A);
            Bp = pointer(B);
            unsafe_copyto!(Bp, Ap, length(A) )
            o = B.^3  + B - B
		end
    end
end


for n in [10000, 1000000, 10000000]
    print("\n n=",n)
    print("\n gpu copy")
    @btime test("gpu copy",$n)
    print("\n gpu")
    @btime test("gpu",$n)
    print("\n cpu")
    @btime test("cpu",$n)
end


#=
N = [1000,10000,100000,1000000,10000000]
gpu_copy = BenchmarkTools.Trial[]
gpu = BenchmarkTools.Trial[]
cpu = BenchmarkTools.Trial[]
for n in reverse(N)
    push!(gpu_copy, @benchmark test("gpu copy",$n))
    push!(gpu, @benchmark test("gpu",$n))
    push!(cpu, @benchmark test("cpu",$n))
end
=#

# t4 = @benchmark test("gpu broadcast")
# t1 = @benchmark test("gpu pin")


