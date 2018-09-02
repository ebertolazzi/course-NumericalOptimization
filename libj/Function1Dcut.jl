#include("../libj/Function1D.jl")


function eval_1Dcut(funct, x0, d, alpha )
  # This method evaluate the function in a point of the "line of cut".
  # $\alpha$ represents the coordinate along the line.
  #print(x0 + alpha * d)
  return evalf(x0 + alpha * d )
end

function eval_D( funct, x0, d, alpha )
  if FD_D
    return Dy(funct, x0, d, alpha )
  else
    return dot( grad( x0 + alpha * d )', d )
  end
end

function eval_DD( funct, x0, d, alpha )
  if FD_DD
    return DDy( funct , x0, d, alpha )
  else
    return dot(hessian( x0 + alpha * d )*d, d )
  end
end

# return the value of the functions computed in x
#function y(funct, x)
#    return eval_1Dcut(x)
#end

# return the value of the first derivative computed in x
function Dy(funct, x0, d, alpha)
    h = (max(1,abs(alpha))*(eps()^(1/3)))
    #print("h = ",h)
    fp = evalf(x0 + (alpha + h) * d )
    #print("fp = ",fp)
    fm = evalf(x0 + (alpha - h) * d )
    #print("fm = ",fm)
    return (fp - fm)/(2*h)
end

function DDy(funct, x0, d, alpha)
    h = max(1 , abs( alpha )) * (eps()^(1/3))
    fp = evalf(x0 + (alpha + h) * d )
    fm = evalf(x0 + (alpha - h) * d )
    fc = evalf(x0 + alpha * d )
    return ( fp + fm -2*fc ) / (h^2)
end
