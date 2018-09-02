include("../libj/FunctionND.jl")
N = 2
M = 1
#exact_solutions = [ 0 ; 0 ]    % one known solution
#guesses         = [ 0.5; 1.0]

function Shubert()
      return x -> x
end

function evalf( x )
      # evaluate function
      check_x(x)
      factor1 = 0.0
      for i = 1 : 5
            factor1 = factor1 + Float64(i) * cos( ( Float64(i) + 1.0 ) * x[1] + Float64(i) )
      end
      factor2 = 0.0
      for i = 1 : 5
            factor2 = factor2 + Float64(i) * cos( ( Float64(i) + 1.0 ) * x[2] + Float64(i) )
      end
      return factor1 * factor2
end

function grad( x )
      # use analitic gradient
      check_x( x )
      g = zeros(Float64, 1, 2 )

      factor1 = 0.0
      df1dx1  = 0.0
      for i = 1 : 5
            factor1 = factor1 + Float64(i) * cos( ( Float64(i) + 1.0 ) * x[1] + Float64(i) )
            df1dx1 = df1dx1 - Float64(i) * ( Float64(i) + 1.0 ) * sin( ( Float64(i) + 1.0 ) * x[1] + Float64(i) )
      end
      factor2 = 0.0
      df2dx2 = 0.0
      for i = 1 : 5
            factor2 = factor2 + Float64(i) * cos( ( Float64(i) + 1.0 ) * x[2] + Float64(i) )
            df2dx2 = df2dx2 - Float64(i) * ( Float64(i) + 1.0 ) * sin( ( Float64(i) + 1.0 ) * x[2] + Float64(i) )
      end

      g[1] = df1dx1 * factor2
      g[2] = factor1 * df2dx2
      return g
end

function hessian( x )
      # use analitic hessian
      check_x( x )
      h = zeros( Float64, 2, 2 )
      factor1 = 0.0
      df1dx1  = 0.0
      df1dx11 = 0.0
      for i = 1 : 5
            factor1 = factor1 + Float64(i) * cos( ( Float64(i) + 1.0 ) * x[1] + Float64(i) )
            df1dx1 = df1dx1 - Float64(i) * ( Float64(i) + 1.0 ) * sin( ( Float64(i) + 1.0 ) * x[1] + Float64(i) )
            df1dx11 = df1dx11 - Float64(i) * ( Float64(i) + 1.0 )^2 * cos( ( Float64(i) + 1.0 ) * x[1] + Float64(i) )
      end
      factor2 = 0.0
      df2dx2 = 0.0
      df2dx22 = 0.0
      for i = 1 : 5
            factor2 = factor2 + Float64(i) * cos( ( Float64(i) + 1.0 ) * x[2] + Float64(i) )
            df2dx2 = df2dx2 - Float64(i) * ( Float64(i) + 1.0 ) * sin( ( Float64(i) + 1.0 ) * x[2] + Float64(i) )
            df2dx22 = df2dx22 - Float64(i) * ( Float64(i) + 1.0 )^2 * cos( ( Float64(i) + 1.0 ) * x[2] + Float64(i) )
      end

      h[1,1] = df1dx11 * factor2
      h[1,2] = df1dx1 * df2dx2
      h[2,1] = df1dx1 * df2dx2
      h[2,2] = factor1 * df2dx22
      return h
end
