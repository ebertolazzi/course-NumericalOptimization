clear all;
close all;
clc;

addpath('../lib');
addpath('../functions');

%r = Rosenbrock();
r = Rastrigin();

disp(r.arity());

%r.contour([-1,3],[-1,2],@(z) log(1+z), 100 );
r.contour([-1,1],[-1,1],@(z) log(1+z), 80)
axis equal ;

search_method   = GoldenSearch();
%search_method   = LinesearchArmijo();
gradient_method = MinimizationGradientMethod(search_method,1e-6,10000,true);

gradient_method.setFunction( r ) ;

%[xs,converged] = gradient_method.minimize( [1;2] ) ;
[xs,converged] = gradient_method.minimize( [0.01;0.0] ) ;

xs
