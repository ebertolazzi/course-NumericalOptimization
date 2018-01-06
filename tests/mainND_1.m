clear all;
close all;
clc;

addpath('../lib');
addpath('../functions');

r = Rosenbrock();
%r = BohachevskyN1();

disp(r.arity());


%ls = 'Armijo' ;
ls = 'Wolfe' ;
%ls = 'GS' ;
switch ls
case 'Armijo'
  linesearch_method = LinesearchArmijo();
case 'Wolfe'
  linesearch_method = LinesearchWolfe();
  linesearch_method.strongWolfe_on();
case 'GS'
  linesearch_method = LinesearchGoldenSection();
  linesearch_method.setMaxIteration( int32(10) ) ;
  linesearch_method.setTolerance(1e-6);
end
  
linesearch_method.debug_on();

%method = 'GRAD' ;
method = 'FR' ;
minimization_method = MinimizationCG( r, linesearch_method ) ;
%minimization_method.selectByName(method);
minimization_method.selectByNumber(14);


minimization_method.setMaxIteration( int32(1000) );
minimization_method.setTolerance(1e-9);
minimization_method.debug_on();
minimization_method.no_FD_D();

fprintf('method = %s\n',minimization_method.activeMethod()) ;

[xs,converged] = minimization_method.minimize( [-1;2] ) ;

subplot(2,1,1) ;
[xmin,ymin,xmax,ymax] = minimization_method.iterRange();
r.contour([xmin,xmax],[ymin,ymax],@(z) z, 80)
axis equal ;
minimization_method.plotIter();

subplot(2,1,2) ;
minimization_method.plotResidual();
xs

fprintf('method = %s\n',minimization_method.activeMethod()) ;

