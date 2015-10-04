"""
Some experiments with Nesterov acceleration on the relaxed basis pursuit:

   min_x 0.5||A*x-b||_2^2 + lambda*||x||_1

See L. Vandenberghe's lecture notes for details.
   http://www.seas.ucla.edu/~vandenbe/236C/
"""

module prox_bp

using PyPlot # only for testing

"""
Basis pursuit via a few proximal gradient methods.
   
   min_x 0.5||A*x-b||_2^2 + lambda*||x||_1

`step=:fixed` - use the keyword arg `step_size` as the fixed step size
`step=:btls` - use a simple backtracking line search, starting with `step_size`
`step=:fista_fixed` - FISTA with a fixed step size `step_size`
`step=:fista_btls1` - FISTA with a simple bactracking line search, starting with `step_size`.
`step=:fista_btls2` - FISTA with a simple bactracking line search, starting with `step_size`.
"""
function bp{T<:Real}(
      A::AbstractArray{T,2},
      b::AbstractArray{T,1};
      lambda::T=1.0,
      step=:fixed,
      step_size::T=1e-3,
      bt_red::T=0.5,
      n_iter::Int=100,
      testing=false)
   
   m,n = size(A)
   tk = step_size # used for *btls 
   x = Array{T}(zeros(n))
   if step in (:fista_fixed, :fista_btls1, :fista_btls2)
      xprev = Array{T}(zeros(n))
      #nu = Array{T}(zeros(n)) # for FISTA reformulation
   end
   
   if testing
      fvals = Array{T}(zeros(n_iter))
   end

   for k in 1:n_iter
      if testing
         fvals[k] = 0.5*norm(A*x-b,2)^2+lambda*norm(x,1)
         #println("f(x) = $(fvals[k])")
      end

      if step == :fixed
         tk = step_size
         grad = A.'*(A*x-b)
         x = prox_l1(x-tk*grad, lambda*tk)

      elseif step == :btls
         x, tk = btls(_x->0.5*norm(A*_x-b,2)^2,
            _x->A.'*(A*_x-b),
            lambda, x, tk, bt_red)

      elseif step == :fista_fixed
         tk = step_size
         y = T(1.0+(k-2.0)/(k+1.0))*x - T((k-2.0)/(k+1.0))*xprev
         xprev[:] = x[:]
         grad = A.'*(A*y-b)
         x = prox_l1(y-tk*grad, lambda*tk)

      elseif step == :fista_btls1
         y = T(1.0+(k-2.0)/(k+1.0))*x - T((k-2.0)/(k+1.0))*xprev
         xprev[:] = x[:]

         x, tk = fista_btls1(_x->0.5*norm(A*_x-b,2)^2,
            _x->A.'*(A*_x-b),
            lambda, y, tk, bt_red)

      else
         error("Unknown step type $(step)")
      end
      
   end

   if testing
      return fvals
   end
end

"""
proximal function for lambda*||x||_1
"""
function prox_l1{T}(x::Array{T,1}, lambda::T)
   lambda > 0 || error("lambda should be positive.")
   p = Array{T}(zeros(size(x)))
   
   inds = x .>= lambda
   p[inds] = x[inds] - lambda
   inds = x .<= -lambda
   p[inds] = x[inds] + lambda

   return p
end


"""
simple backtracking line search
`f` - f(x) = 0.5*||A*x-b||_2^2
`grad` - grad(x) = ∇f(x) = A.'*(A*x-b)
`x` - current x
`t` - step size to start with
`bt_red` - backtracking step size reduction factor
"""
function btls{T<:Real}(
      f::Function,
      grad::Function,
      lambda::T,
      x::Array{T,1},
      t::T,
      bt_red::T)
   
   G(x,t,gradx) = 1/t*(x-prox_l1(x-t*gradx, t*lambda))

   fx = f(x)
   gradx = grad(x)
   tk = t

   Gtx = G(x,tk,gradx)
   upper_bnd = fx - tk*dot(gradx, Gtx) + 0.5*tk*norm(Gtx,2)^2
   xk = x-tk*Gtx
   _f = f(xk)

   while _f > upper_bnd
      tk *= bt_red

      Gtx = G(x,tk,gradx)
      upper_bnd = fx - tk*dot(gradx, Gtx) + 0.5*tk*norm(Gtx,2)^2
      xk = x-tk*Gtx
      _f = f(xk)
   end

   return xk, tk
end


"""
simple backtracking line search method 1 for FISTA
`f` - f(x) = 0.5*||A*x-b||_2^2
`grad` - grad(x) = ∇f(x) = A.'*(A*x-b)
`y` - mix of current x and previous x
`t` - step size to start with
`bt_red` - backtracking step size reduction factor
"""
function fista_btls1{T<:Real}(
      f::Function,
      grad::Function,
      lambda::T,
      y::Array{T,1},
      t::T,
      bt_red::T)
   
   fy = f(y)
   grady = grad(y)

   tk = t
   
   xk = prox_l1(y-tk*grad(y), tk*lambda)
   fxk = f(xk)
   upper_bnd = fy + dot(grady, xk-y) + T(1.0/(2.0*tk))*norm(xk-y,2)^2
   
   while fxk > upper_bnd
      tk *= bt_red

      xk = prox_l1(y-tk*grad(y), tk*lambda)
      fxk = f(xk)
      upper_bnd = fy + dot(grady, xk-y) + T(1.0/(2.0*tk))*norm(xk-y,2)^2
   end 

   return xk, tk
end 




## testing ##
function plot_rel_err(fvals, fs=0, fmt_str="b-")
   
   if fs == 0
      fs = fvals[end]
   end

   rel_err = fvals[1:end-1] - fs
   rel_err /= fs

   semilogy(1:length(rel_err), rel_err, fmt_str)
   xlabel("k")
   ylabel("relative error: (f(xk)-f*)/f*")

end

function run(;timing=false)
   m = 2000; n = 1000;
   #A = 2.0*rand(m,n)-1.0
   A = randn(m,n)
   x = rand(n); b = A*x;
   lambda = 1.0
   fs = 0.5*norm(A*x-b,2)^2+lambda*norm(x,1)
   L = eigmax(A.'*A) # Lipschitz
   
   if timing
      @time fvals_fixed = bp(A, b, lambda=lambda, testing=true, step=:fixed, step_size=1.0/L)
      @time fvals_btls = bp(A, b, lambda=lambda, testing=true, step=:btls, step_size=1e-3)
      @time fvals_fista_fixed = bp(A, b, lambda=lambda, testing=true, step=:fista_fixed, step_size=1.0/L)
      @time fvals_fista_btls1 = bp(A, b, lambda=lambda, testing=true, step=:fista_btls1, step_size=1e-3)
      return
   else
      fvals_fixed = bp(A, b, lambda=lambda, testing=true, step=:fixed, step_size=1.0/L)
      fvals_btls = bp(A, b, lambda=lambda, testing=true, step=:btls, step_size=1e-3)
      fvals_fista_fixed = bp(A, b, lambda=lambda, testing=true, step=:fista_fixed, step_size=1.0/L)
      fvals_fista_btls1 = bp(A, b, lambda=lambda, testing=true, step=:fista_btls1, step_size=1e-3)
   end
  
   clf()
   #fs = minimum(map(x->minimum(x), (fvals_fixed, fvals_btls, fvals_fista_fixed, fvals_fista_btls1)))
   plot_rel_err(fvals_fixed, fs)
   plot_rel_err(fvals_btls, fs, "r-")
   plot_rel_err(fvals_fista_fixed, fs, "g-")
   plot_rel_err(fvals_fista_btls1, fs, "g--")
end

run()
run(timing=true)

end # prox_bp
