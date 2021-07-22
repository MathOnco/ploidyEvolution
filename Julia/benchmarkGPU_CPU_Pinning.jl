using Distributed
using BenchmarkTools
using CUDA
using StatsBase


function test(ttype::String)
    B_ = ones(Float32, 134217728);
    B = cu(B_);
    A = rand(Float32, 134217728);
    A = Mem.pin(Array{eltype(A)}(undef, size(A)));

    Ap = pointer(A);
    Bp = pointer(B);

    idx = sample!(Random.GLOBAL_RNG, 1:134217728, zeros(100))
    idx = Int.(idx)
    idy = sample!(Random.GLOBAL_RNG, 1:134217728, zeros(100))
    idy = Int.(idy)

    for i in 1:20
        if ttype == "cpu"
			o = B_.^8
        elseif ttype == "gpu"
            o = B.^8
        elseif ttype == "gpu broadcast"
            B[1:100] .= 2
            o = B.^8
        elseif ttype == "gpu copy"
            A[idx] = A[idy]
            B = cu(A)
            o = B.^8
        elseif ttype == "gpu pin"
            unsafe_copyto!(Bp, Ap, length(A) )
            o = B.^8
		end
    end
end

t4 = @btime test("gpu broadcast")
t3 = @btime test("cpu")
t2 = @btime test("gpu copy")
t1 = @btime test("gpu pin")
t0 = @btime test("gpu")


