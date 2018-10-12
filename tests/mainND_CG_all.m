clear all;
close all;
clc;

addpath('../lib');
addpath('../functions');

%%
fun_name = 'SchafferF6';
fplot    = @(z) log(1+z);

switch fun_name
case 'Quad'
  r  = Quadratic2D();
  RX = [-1.5 1.5] ;
  RY = [-1.5 1.5] ;
  x0 = [-1; 1 ];
  fplot = @(z) z ;
case 'Rosen'
  r  = Rosenbrock();
  RX = [-1.5 1.5] ;
  RY = [-0.4 1.5] ;
  x0 = [-1.2; 1 ];
case 'Bro'
  r  = Brown_bsf();
  RX = [0,4000] ;
  RY = [-600,600] ;
  x0 = 0.9*[ 10^6  ; 2*10^(-6) ] ;
case 'SchafferF6'
  r = SchafferF6();
  RX = [-20,20] ;
  RY = [-20,20] ;
  x0 = r.guess(int32(1));
  fplot = @(z) z ;
end
%%
%r = Meyer_function;
%disp(r.arity());

ls = 'Wolfe';

switch ls
case 'GS'
  linesearch_method = LinesearchGoldenSection();
  linesearch_method.setTolerance( 1e-9 );
  linesearch_method.setMaxIteration( int32(150) );
case 'Armijo'
  linesearch_method = LinesearchArmijo();
case 'Wolfe'
  linesearch_method = LinesearchWolfe();
end

linesearch_method.debug_on() ;

%min_method = MinimizationGradientMethod(r,linesearch_method);
%min_method = MinimizationBFGS( r, linesearch_method );
min_method = MinimizationCG( r, linesearch_method );

min_method.setMaxIteration( int32(1000) );
min_method.setTolerance(1e-8);
min_method.debug_on();
min_method.no_FD_D();

first = true ;
%for kk=2:14
for kk=2

  min_method.selectByNumber( int32(kk) );
  
  fprintf('\n\nMethod %s\n\n',min_method.activeMethod() );
  [xs,converged] = min_method.minimize( x0 ) ;

  if first
    
    subplot(2,1,1) ;
    r.contour(RX,RY,fplot,80);
    min_method.plotIter();
    axis equal ;

    subplot(2,1,2) ;
    min_method.plotResidual();
    first = false ;
  else
    subplot(2,1,2) ;
    hold on ;
    min_method.plotResidual();
  end
end

xs

fprintf('All done folks!\n\n') ;


