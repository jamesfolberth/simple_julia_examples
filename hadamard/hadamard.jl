
# Code stolen/ported from
# https://github.com/stephenbeckr/SparsifiedKMeans/blob/master/private/hadamard.c


"""
MD doc
"""
function hadamard_jl!{T <: Real}(x::Array{T}, y::Array{T})
# {{{   
   # length must be 2^n
   m = length(x)
   if !ispow2(m) # uses count_ones(m) == 1, which uses popcnt
      error("input length should be a power of 2.")
   end

   if m != length(y)
      error("input length should match output length")
   end

   #XXX Julia's pointer arithmetic does _not_ behave like C's
   #    e.g. xp + 1 moves 1 byte, not sizeof(T) bytes
   #    but we can use unsafe_load/unsafe_store! which will
   #    hide the use of sizeof(T)
   xp = pointer(x,1)
   yp = pointer(y,1)

   k = 0
   for j in 1:2:m-1
      k += 2
      unsafe_store!(yp, unsafe_load(xp, j) + unsafe_load(xp, k), j)
      unsafe_store!(yp, unsafe_load(xp, j) - unsafe_load(xp, k), k)
   end

   bit = 2
   temp = T(0)
   while bit < m
      for j in 1:1:m
         if (bit & (j-1)) == 0
            k = ((j-1) | bit) + 1
            temp = unsafe_load(yp, j)
            unsafe_store!(yp, unsafe_load(yp, j) + unsafe_load(yp, k), j)
            unsafe_store!(yp, temp - unsafe_load(yp, k), k)
         end
      end

      bit <<= 1
   end

# }}}
end


"""
MD doc
"""
function hadamard_c!(x::Array{Float64}, y::Array{Float64})
# {{{
   # This calls Peter Stobbe/Stephen Becker's implementation 
   # in hadamard.so

   # length must be 2^n
   m = length(x)
   if !ispow2(m) # uses count_ones(m) == 1, which uses popcnt
      error("input length should be a power of 2.")
   end

   if m != length(y)
      error("input length should match output length")
   end

   mu = UInt32(m)
   xp = pointer(x) #TODO why use Ref over pointer?
   yp = pointer(y)

   ccall( (:hadamard_apply_vector, "./hadamard.so"), Void,
      (Ref{Float64}, Ref{Float64}, Uint32), yp, xp, mu) # I changed the API :(

# }}}
end


"""
MD doc
"""
function hadamard_c_inplace!(x::Array{Float64}, y::Array{Float64})
# {{{
   # this calls a C version that works in place
   # http://www.musicdsp.org/showone.php?id=18

   # length must be 2^n
   m = length(x)
   if !ispow2(m) # uses count_ones(m) == 1, which uses popcnt
      error("input length should be a power of 2.")
   end

   if m != length(y)
      error("input length should match output length")
   end

   mu = UInt32(m)
   xp = pointer(x) 
   copy!(y,x) # y \gets x
   yp = pointer(y)

   ccall( (:hadamard_apply_vector_inplace, "./hadamard.so"), Void,
      (Ref{Float64}, Uint32), yp, mu)
   
# }}}
end



# Do some timing experiments
############################

# call hadamard_* to compile it ("warm it up")
_x = rand(2^3)
_y = rand(2^3)
hadamard_jl!(_x,_y)
hadamard_c!(_x,_y)
hadamard_c_inplace!(_x,_y)

#x = Array{Float64}([1, 0, 1, 0, 0, 1, 1, 0])
x = rand(2^20)
T = typeof(x[1]) #TODO is there a smarter way to do this?
y = Array{T}(length(x)) # uninitialized array
xh = Array{T}(length(x))


println("hadamard_jl!")
@time hadamard_jl!(x,y)
@time hadamard_jl!(y,xh) # time should be nearly the same
xh /= length(x)

err = norm(x-xh,2)
println("norm(x-xh,2) == $(err)")


println("\nhadamard_c")
y *= 0 # just in case
xh *= 0
@time hadamard_c!(x,y)
@time hadamard_c!(y,xh)
xh /= length(x)

err = norm(x-xh,2)
println("norm(x-xh,2) == $(err)")


println("\nhadamard_c!")
y *= 0 # just in case
xh *= 0
@time hadamard_c_inplace!(x,y)
@time hadamard_c_inplace!(y,xh)
xh /= length(x)

err = norm(x-xh,2)
println("norm(x-xh,2) == $(err)")

