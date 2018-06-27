include("../libj/FunctionND.jl")

N = 4
M = 1
#exact_solutions = [ 1; 1; 1; 1 ];     % one known solution
#guesses         = [ 0.5; 1.0; -0.5; -1.0 ];

function Colville()
      return x -> 100.0 * ( x[2] - x[1]^2 )^2 + ( 1.0 - x[1] )^2 + 90.0 * ( x[4] - x[3]^2 )^2 + ( 1.0 - x[3] )^2 + 10.1 * ( ( x[2] - 1.0 )^2 + ( x[4] - 1.0 )^2 ) + 19.8 * ( x[2] - 1.0 ) * ( x[4] - 1.0 )
end

function evalf( x )
      # evaluate function
      check_x( x )
      return 100.0 * ( x[2] - x[1]^2 )^2 + ( 1.0 - x[1] )^2 + 90.0 * ( x[4] - x[3]^2 )^2 + ( 1.0 - x[3] )^2 + 10.1 * ( ( x[2] - 1.0 )^2 + ( x[4] - 1.0 )^2 ) + 19.8 * ( x[2] - 1.0 ) * ( x[4] - 1.0 )
end

function grad( x )
      # use analitic gradient
      check_x( x )
      g = zeros(Float64, 1, 4 )

      g[1] = 400.0 * x[1]^3 - 400.0 * x[2] * x[1] + 2.0 * x[1] - 2.0
      g[2] = -200.0 * x[1]^2 + 220.2 * x[2] + 19.8 * x[4] - 40.0
      g[3] = -360.0 * x[3] * x[4] + 360.0 * x[3]^3 + 2.0 * x[3] - 2.0
      g[4] = + 180.0 * x[4] - 180.0 * x[3]^2 + 20.2 * x[4] + 19.8 * x[2] - 40.0
      return g
end

function hessian( x )
      # use analitic hessian
      check_x(x)
      h = zeros(Float64, 4, 4 )
      h[1,1] = 1200.0 * x[1]^2 - 400.0 * x[2] + 2.0
      h[1,2] = - 400.0 * x[1]

      h[2,1] = -400.0 * x[1]
      h[2,2] = 220.2
      h[2,4] = 19.8

      h[3,3] = -360.0 * x[4] + 1080.0 * x[3]^2 + 2.0
      h[3,4] = - 360.0 * x[3]

      h[4,2] = 19.8
      h[4,3] = - 360.0 * x[3]
      h[4,4] = 200.2
      return h
end
