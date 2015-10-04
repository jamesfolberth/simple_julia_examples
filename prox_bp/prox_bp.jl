"""
Some experiments with (hybrid) Nesterov acceleration on the basis pursuit:

   min_x 0.5||A*x-b||_2^2 + lambda*||x||_1
"""

module prox_bp

using PyPlot # only for testing

"""
Basis pursuit

`step=:fixed` - use the keyword arg `step_size` as the fixed step size
`step=:btls` - use a simple backtracking line search, starting with `step_size`
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
   x = Array{T}(zeros(n))
   
   if testing
      fvals = Array{T}(zeros(n_iter))
   end

   for k in 1:n_iter
      if testing
         fvals[k] = 0.5*norm(A*x-b,2)^2+lambda*norm(x,1)
         #println("f(x) = $(fvals[k])")
      end

      grad = A.'*(A*x-b)

      if step == :fixed
         tk = step_size
         x = prox_l1(x-tk*grad, lambda*tk)

      elseif step == :btls
         x, tk = btls(_x->0.5*norm(A*_x-b,2)^2,
                  _x->A.'*(A*_x-b),
                  x, step_size, bt_red)

      elseif step == :nest

      else
         error("Unknown step type $(step)")
      end
      
   end

   if testing
      return fvals
   end
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
      x::Array{T,1},
      t::T,
      bt_red::T)
   
   G(x,t,gradx) = 1/t*(x-prox_l1(x-t*gradx, t))

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

function run()
   m = 2000; n = 1000;
   A = randn(m,n);
   b = rand(m);
   lambda = 1.0
   L = eigmax(A.'*A) # Lipschitz
   
   @time fvals_fixed = bp(A, b, lambda=lambda, testing=true, step=:fixed, step_size=1.0/L)
   @time fvals_btls = bp(A, b, lambda=lambda, testing=true, step=:btls, step_size=1e-3)

   clf()
   fs = min(minimum(fvals_fixed), minimum(fvals_btls))
   plot_rel_err(fvals_fixed, fs)
   plot_rel_err(fvals_btls, fs, "r-")
end

run()

end # prox_bp
