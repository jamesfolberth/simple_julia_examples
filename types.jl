"""
Defining a few types
"""
module type_ex

# see julia/base/linalg/factorization.jl and cholesky.jl for an example
abstract Mat{T}
eltype{T}(::Mat{T}) = T

# this should be mutable, otherwise use Diagonal
type DiagMat{T} <: Mat{T}
   diag::Array{T,1} # main diagonal
end

type FullMat{T} <: Mat{T}
   mat::Array{T,2} # full matrix
end

type derived{T, MT<:Mat}
   num::T
   mat::MT
end


# can do some things with metaprogramming
for (t,v) in ((:DiagMat,:diag), (:FullMat,:mat))
   @eval begin
      function op_on_derived{T,MT<:$t}(d::derived{T,MT})
         println(d.mat.$v)
      end
   end
end

# if it is Mat subtype agnostic
function joined_op_on_derived{T,MT<:Mat}(d::derived{T,MT})
   println(d.mat)
end

# specific to DiagMat
function op_on_derived_typed{T,MT<:DiagMat}(d::derived{T,MT})
   println(d.mat)
end


function run()
   
   T = Float64
   #M = DiagMat(Array{T}(rand(3)))
   M = DiagMat(Array{T}([1, 2, 3]))
   #M = FullMat(Array{T}(rand(3,3)))
   #M = FullMat(Array{T}([1 2 3]))
   #M2 = FullMat(Array{T}([1 2 3]))
   #println(typeof(M))
   #println(eltype(M))

   d = derived(T(2),M)
   #d.mat.diag = Array{Int64}([1,2,3,2,3,4])
   #println(d)
   
   #op_on_derived(d)
   #joined_op_on_derived(d)
   op_on_derived_typed(d)

   #println(eltype(d.mat))
end

run()

end # module type_ex
