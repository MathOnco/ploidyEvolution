using BenchmarkTools
using CUDA



function test(ttype::String)
    B_ = ones(Float32, 134217728);
    B = cu(B_);
    A = rand(Float32, 134217728);
    A = Mem.pin(Array{eltype(A)}(undef, size(A)));

    Ap = pointer(A);
    Bp = pointer(B);

    for i in 1:20
        A[:].=1
        if ttype == "cpu"
			o = B_.^8
        elseif ttype == "gpu"
            o = B.^8
        elseif ttype == "gpu index"
            B[1:100] .= 2
            o = B.^8
        elseif ttype == "gpu copy"
            B = cu(A)
            o = B.^8
        elseif ttype == "gpu pin"
            unsafe_copyto!(Bp, Ap, length(A) )
            o = B.^8
		end
    end
end
t4 = @time test("gpu index")
t3 = @time test("cpu")
t2 = @time test("gpu copy")
t1 = @time test("gpu pin")
t0 = @time test("gpu")



