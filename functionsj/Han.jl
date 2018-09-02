include("../libj/FunctionND.jl")
N = 2
M = 1
#exact_solutions = zeros(2,1)
#guesses         = [ 5 ; 3 ]
function Han()
      return x -> x
end

function evalf( x )
      # evaluate function
      check_x(x)
      return x[1]^8 + x[1]^2 + x[1]^2 * x[2]^2 + exp(x[2]^2)
end

function grad( x )
      # use analitic gradient
      check_x(x)
      g   = zeros(Float64, 1, 2 )
      t1  = (x[1]^2)
      t3  = (t1^2)
      t7  = x[2]^2
      t12 = exp(t7)
      g[1] = (8.0 * t3 * t1 * x[1]) + (2.0 * x[1]) + 2.0 * x[1] * t7
      g[2] = (2.0 * t1 * x[2]) + 2.0 * x[2] * t12
      return g
end

function hessian( x )
      # use analitic hessian
      check_x(x)
      h   = zeros(Float64, 2, 2 )
      t1  = x[1]^2
      t2  = t1^2
      t5  = x[2]^2
      t9  = 4.0 * x[1] * x[2]
      t11 = exp(t5)
      h[1,1] = (56.0 * t2 * t1) + 2.0 + 2.0 * t5
      h[1,2] = t9
      h[2,1] = t9
      h[2,2] = (2.0 * t1) + 2.0 * t11 + 4.0 * t5 * t11
      return h
end
