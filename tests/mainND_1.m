clear all;
close all;
clc;

addpath('../lib');
addpath('../functions');

r = Rosenbrock();

disp(r.arity());

%r.contour([-1,3],[-1,2],@(z) log(1+z), 100 );
r.contour([-2,2],[-1,3],@(z) log(1+z), 80)
axis equal ;

search_method = GoldenSearch();
search_method.setMaxIteration( int32(100) ) ;
search_method.setTolerance(1e-10);

%search_method   = LinesearchArmijo();
gradient_method = MinimizationGradientMethod( r, search_method ) ;
gradient_method.setMaxIteration( int32(1000) );
gradient_method.setTolerance(1e-6);
gradient_method.debug_on();

[xs,converged] = gradient_method.minimize( [-1;2] ) ;
gradient_method.plotiter();
xs
