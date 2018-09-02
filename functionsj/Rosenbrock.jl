
# Rosenbrock
# exact_solutions = [1;1];   % one known solution
# guesses         = [-1.2;1] ; % one guess
# workspace()
# p = [-1.1,1.0]
include("../libj/FunctionND.jl")
N = 2
M = 2
function Rosenbrock()
    return x -> [ sqrt(100.0)*(x[2]-x[1]^2); 1.0-x[1] ]
end
function Rosenbrock(x)
    return [ sqrt(100.0)*(x[2]-x[1]^2); 1.0-x[1] ]
end
function evalMap(x)
      # evaluate Rosenbrock (2D) function.
      #X = x[1]
      #Y = x[2]
      return [ sqrt(100.0)*(x[2]-x[1]^2); 1.0-x[1] ]
end

function jacobian( x )
      # use analitic jacobian
      check_x( x )
      #X = x[1]
      return [ -2*sqrt(100.0)*x[1] sqrt(100.0) ; -1.0 0.0 ]
end
function tensor( x )
      # use analitic tensor of second derivative
      T        = zeros(Float64,2,2,2)
      T[1,:,:] = [ -2*sqrt(100.0)  0.0 ; 0.0  0.0 ]
      T[2,:,:] = [ 0.0  0.0 ; 0.0  0.0 ]
      return T
end




# println("evalMap = ",evalMap(p))
# println("evalf = ",evalf(p))
# println("jacobian = ",jacobian(p))
# println("tensor = \n",tensor(p))
# println("grad = ",grad(p))
# println("hessian = ",hessian(p))
# println("FD_grad = ",FD_grad(p))
# println("FD_hessian = ",FD_hessian(p))
# println("FD_jacobian = ",FD_jacobian(p))
# println("FD_tensor = ",FD_tensor(p))


# rewrite of gradient and hessian of the sum of the squares of the map


# Virtual methods definitions
function evalf( x )
      # compute the function as the sum of the squares of the components of the map
      check_x(x)
      F = evalMap( x )
      return 0.5*dot(F,F)
end

function grad( x )
      # gradient of sum f_k(x)^2 --> sum f_k(x) grad_k(x) = 2*F(x)^T J(x)
      check_x(x)
      F = evalMap(x)
      J = jacobian(x)
      return F.' * J
end

function hessian( x )
      # hessian of sum f_k(x)^2 --> 2*(J^T(x) *J(x) + sum f_k(x) Hessian_k(x))
      check_x( x )
      F = evalMap( x )
      J = jacobian( x )
      T = tensor( x )
      H = J.' * J
      for k=1:M
            #H = H + F[k] * squeeze(T[k,:,:])
            H = H + F[k] * T[k,:,:]
      end
      return H
end
