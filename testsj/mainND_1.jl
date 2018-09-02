workspace()

include("../functionsj/Rosenbrock.jl")
include("../libj/Function1Dcut.jl")
include("../libj/MinimizationCG.jl")
r = Rosenbrock()
#p = [1.0;0.0]
#r()
#evalMap([1.0 , 0.0])
LS = "GS"
#LS = "Armijo"
if LS == "GS"
    include("../libj/LinesearchGoldenSection2.jl")
elseif LS == "Armijo"
    include("../libj/LinesearchArmijo2.jl")
else
    print("error: Line Search not defined")
end
#r(p)
#evalf(p)
FD_D = true
FD_DD = true
#eval_1Dcut( r , p , grad( p )' , 0.001 )
#eval_D( r , p , grad( p )' , 0.001 )
#eval_DD( r , p , grad( p )' , 0.001 )
#grad(p)
#step1D( r , p , grad( p )' , 0.001 )
#hessian(p)

#dot(FD_hessian(p)*p, p)
#FD_hessian(p)
method = "GRAD"
#method = "PRP"
#debug_state = true
debug_state = false
#minimization_method = MinimizationCG( r, search_method ) ;
#minimization_method.selectByName(method);
#minimization_method.setMaxIteration( int32(1000) );
#minimization_method.setTolerance(1e-6);
#minimization_method.debug_on();
max_iter = 100000
tol = 1e-12
@time xs,converged = minimize(r ,method, [-0.3;1.5] )
#minimization_method.plotIter();
print(xs)
