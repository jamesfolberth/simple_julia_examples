module module_example

export add_two

function add_two{T <: Real}(x::T)
   x += T(2.)
   return x
end

end
