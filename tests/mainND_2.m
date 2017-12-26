clear all;
close all;
clc;

addpath('../lib');
addpath('../functions');

%r = Rosenbrock();
r = Rastrigin();

disp(r.arity());

r.contour([-1,1],[-1,1],@(z) log(1+z), 80)
axis equal ;


%ls = 'Armijo' ;
%ls = 'Wolfe' ;
ls = 'GS' ;
switch ls
case 'Armijo'
  search_method = LinesearchArmijo();
case 'Wolfe'
  search_method = LinesearchWolfe();
  search_method.strongWolfe_on();
case 'GS'
  search_method = LinesearchGoldenSection();
  search_method.setMaxIteration( int32(10) ) ;
  search_method.setTolerance(1e-6);
end
  
search_method.debug_on();

method = 'gradient' ;
%method = 'PR' ;
switch method 
case 'gradient'
  minimization_method = MinimizationGradientMethod( r, search_method ) ;
case 'PR'
  minimization_method = MinimizationPolackRibiere( r, search_method ) ;
end

minimization_method.setMaxIteration( int32(1000) );
minimization_method.setTolerance(1e-6);
minimization_method.debug_on();

[xs,converged] = minimization_method.minimize( [0.5;0.8] ) ;
minimization_method.plotiter();
xs
