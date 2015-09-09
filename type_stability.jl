
# example from
# https://en.wikibooks.org/wiki/Introducing_Julia/Types#Type_stability
#function ts1(n)
#   s = 0 # this is an instance of an Intger (concrete) subtype
#   for i in 1:n
#      s += s/i # i and s are both Integers (Int64 on my machine),
#               # but this gets promoted to Float64 each iteration.
#               # this is _not_ "type stable"
#               # Julia doesn't know that 0 divided by (almost) anything is 
#               # zero, so it can't optimize out the code.
#   end
#end
#
#function ts2(n)
#   s = 0.0 # this is an instance of an AbstractFloat subtype
#   for i in 1:n
#      s += s/i # Julia can make pure primitive code for this!
#               # this code _is_ "type stable" and gives a huge performance
#               # gain
#   end
#end

# Okay, so the above example isn't exactly fair.  Since Julia can see
# that s is derived from Float and is equal to 0.0, it knows that s/i == 0.0
# no matter the non-zero value of i.  It actually compiles out floating point
# division, hence the enormous speedup.  In fact, I think it compiles out the
# whole loop!  Check the LLVM bitcode and see for yourself.

# Perhaps a more fair example?
function ts1(n)
   s = 1 # this is an instance of an Intger (concrete) subtype
   for i in 1:n
      s += s/i # i and s are both Integers (Int64 on my machine),
               # but this gets promoted to Float64 each iteration.
               # this is _not_ "type stable"
   end
   return s
end

function ts2(n)
   s = 1.0 # this is an instance of an AbstractFloat subtype
   for i in 1:n
      s += s/i # Julia can make pure primitive code for this!
               # this code _is_ "type stable" and gives a good performance 
               # improvement.
               # However, since s is initialized to 1.0, we can't compile
               # away all the floating point divisions.  This is probably a
               # better comparison.
   end
   return s
end

@time ts1(10000000)
@time ts2(10000000)

# check out the LLVM bitcodes for both methods
# the bitcode for ts2 is a whole lot cleaner than for ts1, since
# Julia cannot make good assumptions about the variables in ts1
#@code_llvm ts1(100)
#@code_llvm ts2(100)

