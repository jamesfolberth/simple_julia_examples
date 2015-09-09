
# will this by type stable?
function fT{T <: Number}(a::T, n::Integer)
   s = T(0) # yes, because of this conversion
   for i in 1:n
      s += a
   end
   return s
end

function fFloat64(a::Float64, n::Int64)
   s = Float64(0)
   for i in 1:n
      s += a
   end
   return s
end

function f_nonTS{T <: Number}(a::T, n::Integer)
   s = 0 # this makes us not type stable
   for i in 1:n
      s += a
   end
   return s
end

@time fT(1, 10000000)
@time fFloat64(1., 10000000)

@time f_nonTS(1, 10000000)
@time f_nonTS(1., 10000000)

