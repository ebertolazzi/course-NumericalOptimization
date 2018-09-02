
function minimize(f, interval )
  a = interval[1]
  b = interval[2]
  tau      = (sqrt(5)-1)/2
  tol      = 1e-5
  max_iter = 100
  dlen   = tau*(b-a)
  lambda = b-dlen
  mu     = a+dlen
  fl     = f(lambda);
  fm     = f(mu);
  tolerance = tol*(b-a)
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
      fm     = f(mu)
    else
      # select interval [a,mu]
      b      = mu
      mu     = lambda
      fm     = fl
      dlen   = tau*(b-a)
      lambda = b-dlen
      fl     = f(lambda)
    end
  end
  return a,b
end
