using Distributed
using BenchmarkTools
using CUDA
using StatsBase
using Random


function test(ttype::String)
    B_ = ones(Float32, 1034200);
    B = cu(B_);
    A = rand(Float32, length(B_));
    A = Mem.pin(Array{eltype(A)}(undef, size(A)));

    Ap = pointer(A);
    Bp = pointer(B);

    idx = sample!(Random.GLOBAL_RNG, 1:length(B_), zeros(100))
    idx = Int.(idx)
    idy = sample!(Random.GLOBAL_RNG, 1:length(B_), zeros(100))
    idy = Int.(idy)

    for i in 1:20
        if ttype == "cpu"
			o = B_.*B_ + B_ - B_
        elseif ttype == "gpu"
            o = B.*B + B - B
        elseif ttype == "gpu broadcast"
            B[1:100] .= 2
            o = B.*B + B - B
        elseif ttype == "gpu copy"
            #A[idx] = A[idy]
            B = cu(A)
            o = B.*B + B - B
        elseif ttype == "gpu pin"
            unsafe_copyto!(Bp, Ap, length(A) )
            o = B.*B + B - B
		end
    end
end

t2 = @btime test("gpu copy")
t0 = @btime test("gpu")
t3 = @btime test("cpu")
t4 = @btime test("gpu broadcast")
t1 = @btime test("gpu pin")


