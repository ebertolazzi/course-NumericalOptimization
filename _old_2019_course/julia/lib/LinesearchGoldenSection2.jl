include("ForwardBackward2.jl")

function minimize(funct, x0, d, interval )
  a = interval[1]
  b = interval[2]
  tau      = (sqrt(5)-1)/2
  tol_ls      = 1e-4
  max_iter = 100
  dlen   = tau*(b-a)
  lambda = b-dlen
  mu     = a+dlen
  fl     = eval_1Dcut( funct, x0, d, lambda )
  fm     = eval_1Dcut( funct, x0, d, mu )
  tolerance = tol_ls*(b-a)
  for k = 1:max_iter
    if (b-a) < tolerance
      # print("find a solution\n\n")
      break
    end
    if fl > fm
      # select interval [lambda,b]
      a      = lambda
      lambda = mu
      fl     = fm
      dlen   = tau*(b-a)
      mu     = a+dlen
      fm     = eval_1Dcut( funct, x0, d, mu )
    else
      # select interval [a,mu]
      b      = mu
      mu     = lambda
      fm     = fl
      dlen   = tau*(b-a)
      lambda = b-dlen
      fl     = eval_1Dcut( funct, x0, d, lambda)
    end
  end
  return a,b
end


function search( funct, x0, d, alpha_guess )
  aLO,aHI,fLO,fHI,ierr = ForwardBackward(funct, x0, d, alpha_guess )
  ok = ierr >= 0
  if ok
    a,b = minimize(funct, x0, d, [aLO, aHI] )
    alpha = (a+b)/2.0
    return alpha,ok
  else
    alpha = alpha_guess/100.0
    return alpha,ok
  end
end
