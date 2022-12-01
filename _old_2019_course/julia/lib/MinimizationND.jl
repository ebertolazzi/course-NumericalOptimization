#include("../libj/LinesearchGoldenSection2.jl")


function step1D( funct, x0, d, alpha_guess )

      if norm(d,Inf) == 0
            println("MinimizationND, bad direction d == 0")
      end
      # search an interval for minimization
      alpha,ok = search( funct, x0, d, alpha_guess )

      # check error
      if !ok
            println("MinimizationGradientMethod:step1D, linesearch failed")
      end
      # advance
      x1 = x0 + alpha * d
      return x1, alpha
end
