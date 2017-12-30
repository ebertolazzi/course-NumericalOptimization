clear all;
close all;
clc;

addpath('../lib');
addpath('../functions');

fun_name = 'Quad';

switch fun_name
    case 'Quad'
    r = Quadratic2D();
    disp(r.arity());
    r.contour([-1.5 1.5],[-1.5 1.5],@(z) log(1+z), 80)
    axis equal ;
    x0 = [-1; 1 ];
  case 'Rosen'
    r = Rosenbrock();
    disp(r.arity());
    r.contour([-1.5 1.5],[-0.4 1.5],@(z) log(1+z), 80)
    axis equal ;
    x0 = [-1; 1 ];
  case 'Bro'
    r = Brown_bsf();
    disp(r.arity());
    r.contour([0,4000],[-600,600],@(z) log(1+z), 80)
    axis equal ;
    x0 = 0.9*[ 10^6  ; 2*10^(-6) ];
end


search_method   = GoldenSearch();
%search_method = LinesearchArmijo();

search_method.debug_on();

%dir_method = MinimizationGradientMethod(r,search_method);
dir_method = MinimizationBFGS(r,search_method);

dir_method.setMaxIteration( int32(1000) );
dir_method.setTolerance(1e-6);
dir_method.debug_on();

[xs,converged] = dir_method.minimize( x0 ) ;
dir_method.plotiter();

xs


