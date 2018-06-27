include("../libj/FunctionND.jl")
N = 2
M = 1

function Leon()
    return x -> 100.0*(x[2] - x[1] * x[1] * x[1])^2 + (1.0 - x[1])^2
end

function evalf( x )
      # evaluate function
      check_x(x)
      f1 = x[2] - x[1] * x[1] * x[1]
      f2 = 1.0 - x[1]
      return (100.0 * f1 * f1 + f2 * f2)
end

function grad( x )
      # use analitic gradient
      check_x(x)
      g    = zeros( Float64, 1, 2 )
      g[1] = - 600.0 * ( x[2] - x[1]^3 ) * x[1] * x[1] - 2.0 * ( 1.0 - x[1] )
      g[2] = 200.0 * ( x[2] - x[1]^3 )
      return g
end

function hessian( x )
      # use analitic hessian
      check_x(x)
      h = zeros( Float64, 2, 2 )
      h[1,1] = - 1200.0 * x[1] * x[2] + 3000.0 * x[1]^4 + 2.0
      h[1,2] = - 600.0 * x[1] * x[1]
      h[2,1] = - 600.0 * x[1] * x[1]
      h[2,2] = 200.0
      return h
end
