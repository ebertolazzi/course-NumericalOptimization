include("../libj/MinimizationND.jl")

function PRP( g1, g0, d0 )
  return dot(g1,g1-g0)/dot(g0,g0)
end

function betaEval( method, g1, g0, d0 )
  # Compute the value of `beta` for the selected nonlinear CG method
  if method == "GRAD"  # standrd gradient method
    beta = 0.0
  elseif method == "PRP"
    beta = PRP( g1, g0, d0 )
  else
    print("MinimizationCG, betaEval no method selected (this should not be happen)")
  end
end

    # function self = MinimizationCG( fun, ls )
    #   %
    #   % fun = function to be minimized
    #   % ls  = linesearch used, can be LinesearchArmijo, LinesearchWolfe, LinesearchGoldenSection
    #   self@MinimizationND( fun, ls );
    #   self.method_names    = {'GRAD','HS','HS+','FR','PRP','PRP+','CD','LS','DY','N','HyTAS','HyHS','HyNG','HyDY'};
    #   self.method          = self.method_names{1};
    #   self.direction_short = 1e-3;
    #   self.angle_too_small = pi/180; % 1 degree
    # end
direction_short = 1e-3
angle_too_small = pi/180




function  minimize(funct, method, x0 )
    # generic onjugate gradient minimization algorithm
    # x0 = initial point

    xs    = copy(x0)
    alpha = 1.0
    converged = false
    g0   = copy(x0)
    d0   = copy(x0)
    for iter=1:max_iter
        # gradient of the function
        g1 = grad( xs ).'
        # check if converged
        nrm_g1 = norm(g1,Inf)
        converged = nrm_g1 < tol
        if converged
            if debug_state
                println("solution found, ||grad f||_inf = ",nrm_g1," <", tol )
            end
            println("\nIter : ",iter,"\n\nSolution found, in : ",xs,"\n\n||grad f||_inf = ",nrm_g1," <", tol,"\n\n" )
            return xs, converged
            #break
        end

        # only for debug
        if debug_state
            println("iter = ",iter," ||grad f||_inf = ", nrm_g1 )
        end

        # build search direction

        if iter == 1
            #print("pass here!")
            #g0 = 0
            d1   = -g1 # first iteration, search direction is -gradient
            beta = 0.0
        else
            beta = betaEval(method, g1, g0, d0 )
            d1   = -g1 + beta*d0
            # check if the direction is descending, if direction is >= 89 degree from gradient reset to gradient directions
            # use >= to catch d1 == 0
            nrm_d1 = norm(d1)
            if dot( d1, g1 ) >= -nrm_d1 * nrm_g1 * (pi-angle_too_small)
                if debug_state
                    println("direction angle about 90 degree, reset direction search!")
                end
                d1 = -g1 # reset direction
            elseif nrm_d1 <= direction_short * nrm_g1
                if debug_state
                    println("direction length too short, reset direction search!" )
                end
                d1 = -g1 # reset direction
            end
        end
        if norm(d1,Inf) == 0
            println("MinimizationCG, bad direction d == 0")
        end

        # minimize along search direction
        xs,alpha = step1D(funct, xs, d1, alpha )

        if debug_state
            println(" alpha = ",alpha,", beta = ", beta)
        end
        # save old gradient and direction
        g0 = copy(g1)
        d0 = copy(d1)
    end
    return xs, converged
end
