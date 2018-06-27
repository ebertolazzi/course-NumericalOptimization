include("../libj/FunctionND.jl")
N = 2
M = 1

#exact_solutions = [ 3.0; 2.0]     % one known solution
#guesses         = [ -1.3; 2.7]

function HB()
      return x ->  (x[1]^2 + x[2] - 11.0 )^2 + ( x[1] + x[2]^2 - 7.0 )^2
end

function evalf( x )
      # evaluate function
      check_x( x )
      return (x[1]^2 + x[2] - 11.0 )^2 + ( x[1] + x[2]^2 - 7.0 )^2
end

function grad( x )
      # use analitic gradient
      check_x( x )
      g = zeros(Float64, 1, 2 )
      g[1] = 4.0 * ( x[1]^2 + x[2] - 11.0 ) * x[1] + 2.0 * ( x[1] + x[2]^2 - 7.0 )
      g[2] = 2.0 * ( x[1]^2 + x[2] - 11.0 ) + 4.0 * ( x[1] + x[2]^2 - 7.0 ) * x[2]
      return g
end

function hessian( x )
      # use analitic hessian
      check_x( x )
      h = zeros(Float64, 2, 2 )
      h[1,1] = 8.0 * x[1] * x[1] + 4.0 * ( x[1] * x[1] + x[2] - 11.0 ) + 2.0
      h[1,2] = 4.0 * x[1] + 4.0 * x[2]
      h[2,1] = 4.0 * x[1] + 4.0 * x[2]
      h[2,2] = 2.0 + 8.0 * x[2] * x[2] + 4.0 * ( x[1] + x[2] * x[2] - 7.0 )
      return h
end
