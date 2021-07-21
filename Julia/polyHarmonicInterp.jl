module PolyharmonicInterpolation
using LinearAlgebra

export PolyharmonicInterpolator,polyharmonicSpline

struct PolyharmonicInterpolator
	knots::AbstractArray
	values::Vector{Float64}
	nPts::Int
	dim::Int
	order::Int
	w::Vector{Float64}
	v::Vector{Float64}
end

PolyharmonicInterpolator(knots::AbstractArray{T},values::Vector{T},
order::Int=2) where T <: Real = begin

	# Check to make sure knots is N x d and values is N x 1
	if size(knots)[1] != length(values)
		error("knots and values length must agree")
	end

	if length(size(knots)) == 1
		nPts = length(knots)
		dim = 1
	else
		nPts,dim = size(knots)
	end

	# initialize w and v
	w = Vector{Float64}(undef,nPts)
	v = Vector{Float64}(undef,dim+1)

	# Fill npts and dim using the dimension of knots
	computeWeights(PolyharmonicInterpolator(knots,values,nPts,dim,order,w,v))

end

function polyharmonicRBF(r::Float64,k::Int)


	if r >= 1
		val = (k % 2 == 0) ? r^k*log(r) : r^k
	else
		val = (k % 2 == 0) ? r^(k-1)*log(r^r) : r^k
	end

	return val

end

function createBlockMatrix(p::PolyharmonicInterpolator)

	# Create the matrix of differences between the knots
	dist = [norm(p.knots[i,:] - p.knots[j,:]) for i in 1:p.nPts, j in 1:p.nPts]

	# Create the matrix A
	A = map(r->polyharmonicRBF(r,p.order),dist)

	# Create the matrix B
	B = [ones(p.nPts) p.knots];

	return A,B


end

function computeWeights(p::PolyharmonicInterpolator)

	# Create block matrices used to solve for weights
	A,B = createBlockMatrix(p)

	# Make block matrix for the symmetric linear system
	blockMat = [A B; transpose(B) zeros(p.dim+1,p.dim+1)]

	# Get weights nPts and number of dimensions
	weights = blockMat \ [p.values; zeros(p.dim+1,1)]

	p.w .= weights[1:p.nPts]
	p.v .= weights[p.nPts+1:end]

	return p

end

function polyharmonicSpline(p::PolyharmonicInterpolator,interpPoints::AbstractArray)

	# Check to see the dimensions of the interpPoints array agree
	let
		nInterpDims = nothing
		try
			_, nInterpDims = size(interpPoints)
		catch
			nInterpDims = 1
		end

		if p.dim != nInterpDims
			error("dimension of knots and points to be interpolated must agree")
		end
	end

	# Initialize array
	interpValues = zeros(size(interpPoints)[1])

	for (i,point) in enumerate(eachrow(interpPoints))

		tmppoint =[1.0;point]

		# Get distance to all scattered points
		dist = sqrt.(sum((p.knots .- point').^2,dims=2))

		# reduce to normal size (Julia returns a 2 dimensional row array for some reason...)
		dist = reduce(vcat,dist)

		# Get distance between basis and point of interest and compute the
		# polyharmonic basis
		phi = map(r->polyharmonicRBF(r,p.order),dist)

		interpValues[i] = p.w'phi + p.v'tmppoint

	end

	return interpValues

end

end
