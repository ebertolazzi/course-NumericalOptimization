
function check_x( x )
      szdim = length(size(x))
      #print(size(x))
      if szdim != 2
            print("FunctionND, number of dimension of x = ",szdim," expected 2\n")
      end
      n = size(x,1)
      m = size(x,2)
      if !( n == N && m == 1 )
            print("FunctionND, size(x) = ",n," x ",m,", expected ",N," x 1\n")
      end
end


function FD_grad( x )
      # finite difference approximation of the gradient
      h  = max.(1,abs.(x))*eps()^(1/3)
      xp = copy(x)
      xm = copy(x)
      #@printf "xp[1] = %.10f\n" xp[1]
      #@printf "xp[2] = %.10f\n" xp[2]
      #@printf "xm[1] = %.10f\n" xm[1]
      #@printf "xm[2] = %.10f\n" xm[2]
      g  = zeros(Float64,1,N)
      for k=1:N
            #println("h = ",h)
            xp[k] = x[k]+h[k]

            #@printf "xp = %.10f\n" xp[k]
            xm[k] = x[k]-h[k]
            #@printf "xm = %.10f\n" xm[k]
            fp    = evalf(xp)
            #println("fp = ",typeof(fp))
            fm    = evalf(xm)
            #println("fm = ",fm)
            g[k]  = (fp-fm)./(2*h[k])
            #println("g = ",(xm-xp))
            xp[k] = x[k]
            xm[k] = x[k]
      end
      return g
end

function FD_hessian( x )
      # finite difference approximation of the hessian
      # Baseed on a code by Brendan C. Wood
      # Copyright (c) 2011, Brendan C. Wood <b.wood@unb.ca>
      H = zeros(Float64,N,N)
      h = max.(1,abs.(x))*eps()^(1/3) #%0.1*ones(3,1); % ricetta di cucina, attenzione x e' un vettore
      # for each dimension of objective function
      for i=1:N
            # derivative at first point (left)
            x1    = copy(x)
            x1[i] = x[i] - h[i]
            df1   = grad(x1)

            # derivative at second point (right)
            x2    = copy(x)
            x2[i] = x[i] + h[i]
            df2   = grad(x2)

            # differentiate between the two derivatives
            d2f = (df2-df1) ./ (2*h[i])

            # assign as column i of Hessian
            H[:,i] = d2f
      end

      # Make H symmetric numerically, this is not mandatory but could help
      return 0.5*(H+H.')
end

function FD_jacobian( x )
      # finite difference approximation of the jacobian
      check_x(x)
      h  = max.(1,abs.(x))*eps()^(1/3)
      xp = copy(x)
      xm = copy(x)
      J  = zeros(Float64,M,N)
      for k=1:N
            xp[k]  = x[k]+h[k]
            xm[k]  = x[k]-h[k]
            #fp     = evalMap(xp)
            #fm     = evalMap(xm)
            fp     = evalf(xp)
            fm     = evalf(xm)
            J[:,k] = (fp-fm)./(2*h[k])
            xp[k]  = x[k]
            xm[k]  = x[k]
      end
      return J
    end

function FD_tensor( x )
      # finite difference approximation of the second derivative
      # T(k,:,:) is the hessian of the k-th components of the map
      check_x(x)
      T = zeros(Float64,M,N,N)
      h = max.(1,abs.(x))*eps()^(1/3) ##0.1*ones(3,1); % ricetta di cucina, attenzione x e' un vettore
      # for each dimension of objective function



      for i=1:N
            # derivative at first point (left)
            x1    = copy(x)
            x1[i] = x[i] - h[i]
            jf1   = jacobian(x1)

            # derivative at second point (right)
            x2    = copy(x)
            x2[i] = x[i] + h[i]
            jf2   = jacobian(x2)

            # differentiate between the two derivatives
            j2f = (jf2-jf1) ./ (2*h[i]) #% M x N

            # assign as column i of Hessian
            T[:,i,:] = j2f
      end

      # Make T symmetric numerically, this is not mandatory but could help
      for k=1:M
            T[k,:,:] = 0.5*(T[k,:,:]+T[k,:,:].')
      end
      return T
end
