
function ForwardBackward( f, alpha_guess )
      # find alpha_min <= alpha0 < alpha1 <= alpha_max such that
      # alpha0 satify Armijo and alpha1 DO NOT satisfy Armijo
      # ierr = 0  interval found
      # ierr = 1  both alpha0, alpha1 satisfy Armijo f(alpha0) <  f(alpha1)
      # ierr = 2  both alpha0, alpha1 satisfy Armijo f(alpha0) >= f(alpha1)
      # ierr = -1 both alpha0, alpha1 DO NOTY satisfy Armijo
      # ierr = -2 Df0 >= 0
      c1 = 0.1
      c2 = 0.5
      tau_LS = 1.5
      tau_acc = 1.2
      alpha_min = 1e-10
      alpha_max = 1e10
      dumpMin = 0.05
      dumpMax = 0.951
      alpha_epsi = eps()^(1/3)
      debug_status = false
      #f = funct
      # correct alpha_guess into the required interval
      alpha0 = max(alpha_min,min(alpha_max,alpha_guess))
      # compute initial value and derivative
      #f        = @(a) self.fun1D.eval(a) ;
      #df       = @(a) self.fun1D.eval_D(a) ;
      f0  = y(f, 0)
      Df0 = Dy(f, 0)
      #println("f(0) = ", f0)
      #println("Df(0) = ", Df0)
      # if not decreasing issue an error (if in debug also plot the function)
      #if Df0 >= 0
      #  println("Error Df0 = ", Df0, " must be negative")
      #end
      # initialize search parameters
      c1Df0 = c1*Df0
      tauf  = tau_LS
      ierr  = 0
      # decide if do forward or backward search
      fa0 = y(f, alpha0)
      if fa0-f0 <= alpha0 * c1Df0
        # satisfy Armijo --> forward search

        # if increasing minima in [0,alpha0]
        if Dy(f, alpha0) >= 0
          alpha1 = alpha0
          fa1 = fa0
          alpha0 = 0
          fa0 = f0
          ierr = 1
          return alpha0,alpha1,fa0,fa1,ierr
        end

        # not increasing search next interval
        alpha1 = tauf * alpha0
        fa1 = y(f, alpha1)
        # continue to loop until satify Armijo and is non increasing
        while true
          # check if found point that violates Armijo
          if fa1-f0 > alpha1 * c1Df0
                 return alpha0,alpha1,fa0,fa1,ierr
          end
          # if increasing break
          if fa1 > fa0
                ierr = 1
                return alpha0,alpha1,fa0,fa1,ierr
          end
          if Dy(f, alpha1) >= 0
                ierr = 2
                return alpha0,alpha1,fa0,fa1,ierr
          end

          # check if interval become too large
          if tauf * alpha0 >= alpha_max
            # last interval to check
            alpha0 = alpha1
            fa0 = fa1
            alpha1 = alpha_max
            fa1 = y(f, alpha1)
            if fa1 > f0 + alpha1 * c1Df0
                  return alpha0,alpha1,fa0,fa1,ierr
            end # found!
            if fa1 > fa0
                  ierr = 1
            else
                  ierr = 2
            end
            return alpha0,alpha1,fa0,fa1,ierr
          end
          # prepare for next interval
          alpha0 = alpha1
          fa0 = fa1
          alpha1 = tauf * alpha0
          fa1 = y(f, alpha1)
          tauf = tauf * tau_acc # update tau factor
        end
        # never pass here
      else
        # DO NOT satisfy Armijo --> backward search
        alpha1 = alpha_guess
        fa1 = fa0
        alpha0 = alpha1/tauf
        fa0 = y(f, alpha0)
        while fa0-f0 > alpha0 * c1Df0 # DO NOT satisfy Armijo
          alpha1 = alpha0
          fa1 = fa0
          if alpha1 <= alpha_min * tauf
            # force termination, last check
            alpha0 = alpha_min
            fa0 = y(f, alpha0)
            if fa0-f0 > alpha0 * c1Df0
                  ierr = -1
            end
            return alpha0,alpha1,fa0,fa1,ierr
          end
          alpha0 = alpha1/tauf
          fa0 = y(f, alpha0)
          tauf = tauf * tau_acc # update tau factor
        end
        # exiting from loop satisfying Armijo
        #return alpha0,alpha1,fa0,fa1,ierr
      end
end
