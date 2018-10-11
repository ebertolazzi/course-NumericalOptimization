clear all;
close all;
clc;

addpath('../lib');
addpath('../functions');

%%
fun_name = 'Boha2';

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
  x0 = r.guesses();
case 'Brown_badly_scaled'
  r = Brown_bsf();
  disp(r.arity());
  r.contour( [10^6 -2000,10^6 + 2000],[-1000,1000],@(z) log(1+z), 80)
  %axis equal ;
  x0 = r.guesses();
case 'SchafferF6'
  r = SchafferF6();
  disp(r.arity());
  r.contour([-20,20],[-20,20],@(z) z, 80)
  axis equal ;
  x0 = r.guesses();
case 'Boha1'
  r = BohachevskyN1();
  disp(r.arity());
  r.contour([-4 4],[-4 4],@(z) log(1+z), 80)
  axis equal ;
  x0 = r.guesses();
case 'Boha2'
  r = BohachevskyN2();
  disp(r.arity());
  r.contour([-4 4],[-4 4],@(z) log(1+z), 80)
  axis equal ;
  x0 = r.guesses();
case 'Boha3'
  r = BohachevskyN3();
  disp(r.arity());
  r.contour([-4 4],[-4 4],@(z) log(1+z), 80)
  axis equal ;
  x0 = r.guesses();
case 'Hilb'
  r = Hilbert();
  disp(r.arity());
  r.contour([-4 4],[-4 4],@(z) log(1+z), 80)
  axis equal ;
  x0 = r.guesses();
end

%%
clear r x0
clf
r = Brown_Dennis_function(20);
%r = Powell_singular_function;
%r = Rosenbrock();
x0 = r.guesses;

%%

search_method  = LinesearchGoldenSection();
%search_method = LinesearchArmijo();
%search_method = LinesearchWolfe();

search_method.debug_on();

%dir_method = MinimizationGradientMethod(r,search_method);
%dir_method  = MinimizationBFGS(r,search_method);
%dir_method = MinimizationQN(r,search_method);
dir_method = MinimizationCG(r,search_method);

dir_method.setMaxIteration( int32(1000) );
dir_method.setTolerance(1e-6);
dir_method.debug_on();
%dir_method.selectByName('BFGS');
%dir_method.selectByName('DFP');
%dir_method.selectByName('PSB');

[xs,converged] = dir_method.minimize( x0 ) ;
dir_method.plotIter();

fprintf('\nThe solution found by the algorithm is:\n');
fprintf('       %d\n',xs)
exact = r.exact_solutions;
fprintf('\nThe analytical solution is:\n');
fprintf('       %d\n',exact)



