"""
Some simple tests to help determine if your gradient computations are sensible.

   Tests: taylor_test
"""

module grad_test

# take a look at the doc for cholfact for an example of documentation
"""
`taylor_test(f, g, x0, [;pert=1e-3, n_iter=5])`

Try a few simple Taylor series based tests.  See Stephen's gradientTest.pdf for
the derivation.

The function `f` should map R^n -> R^k and have syntax `f(x)`.  To use a
function with additional parameters, use an anonymouse function: 
`x -> f(x, args...)`.  The function `g` should compute the gradient of `f`.  
`x0` is a point in the domain of `f` and `g`.  `pert` controls how large of a
perturbation from `x0` we use when computing Taylor series terms.
"""
function taylor_test{T<:Real}(
      f::Function,
      g::Function,
      x0::DenseArray{T};
      pert=T(1e-3),
      n_iter=5)
   
   # helper
   my_norm(x::T) = abs(x)
   my_norm(x::DenseArray{T}) = sqrt((x[:]'*x[:])[1])
   
   # build a perturbation to x0
   ep = Array{T}(pert*rand(size(x0)))

   # we'll scale down the perturbation during the iterations below
   const C = T(10)
   ep *= C^2 
   println("Running gradient test:")
   println("||x1-x0||,    little o,     gradient IP,  taylor error") 
   
   f0 = f(x0); g0 = g(x0)

   # check that things are real
   if !(eltype(f0) <: Real) || !(eltype(g0) <: Real)
      warn("This hasn't been tested for complex functions!")
   end
   
   for i in 1:n_iter
      ep /= C
      x1 = x0 + ep
      f1 = f(x1); g1 = g(x1)
      
      if my_norm(f1-f0) < eps(f0)*T(1e4)
         println("Difference ||f(x1)-f(x0)|| = $(my_norm(f1-f0)) too small.")
         continue
      end
      
      # Is the Talor truncation valid? This term should be small, but 
      # it wont necessarily go to zero with ep.
      little_o_small = maximum(abs(my_norm(ep)^2 ./ (T(2)*(f1-f0))))

      # the (g0[:]+... will compute the usual inner product on 
      #  C^n: <x,y> = x'*y
      #  S(C)^n: <X,Y> = tr(X'*Y) = vec(X)'*vec(Y)
      #
      # This is the term involving the IP of gradients derived in
      # Stephen's gradientTest.pdf.  This term should go to zero.
      grad_IP = ((g0[:]+g1[:])'*ep[:]) ./ (T(2)*(f1-f0))
      grad_IP_small = maximum(abs(grad_IP-T(1)))

      # normalized first-order taylor error
      taylor_er = maximum(abs(f0 + g0[:]'*ep[:] - f1) / my_norm(ep))
 
      for term in 
         (my_norm(ep), little_o_small, grad_IP_small, taylor_er)
         print(@sprintf("%-14.2E", term))
      end
      println("")
   end
end


"""
`simple_taylor_test()`

Test `taylor_test` on

    f(x) = 0.5||A*x-b||_2^2
    g(x) = ∇f(x) = A'*(A*x-b)
"""
function simple_taylor_test()
   const m = 4
   const n = 5
   const A = rand(m,n)
   const b = rand(m,1)
   function f1(x, A, b)
      return 0.5*norm(A*x-b,2)^2
   end
   f1(x) = f1(x, A, b)
   
   function g1(x, A, b)
      return A.'*(A*x-b)
   end
   g1(x) = g1(x, A, b)
   
   x = rand(n)
   taylor_test(f1, g1, x)
end

end
