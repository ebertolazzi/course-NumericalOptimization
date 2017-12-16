clear all;
addpath('../lib');

% define the object function fs, fc
fs = SinFun() ;
fc = CosFun() ;

x = 0:2*pi/1000:2*pi ;

subplot( 2, 1, 1);
hold off ;
disp('compute with fs') ;
plot( x, fs.eval(x), ...
      x, fs.eval_D(x), ...
      x, fs.eval_DD(x)) ;
title('function fs') ;
hold on ;

subplot( 2, 1, 2);
hold off ;
disp('compute with fc') ;
plot( x, fc.eval(x), ...
      x, fc.eval_D(x), ...
      x, fc.eval_DD(x)) ;
title('function fc') ;
hold on ;

disp('do golden search') ;
% initialize Golden Search object m
m = GoldenSearch();
% attach function fs to the miminization object
m.setFunction(fs);

[a,b] = m.minimize(pi,2*pi);
fa = fs.eval(a) ;
fb = fs.eval(b) ;

subplot( 2, 1, 1);
plot( x, fs.eval(x) ) ;
plot( pi, fs.eval(pi),     'ok' ) ;
plot( 2*pi, fs.eval(2*pi), '*k' ) ;
plot( a, fa, 'or' ) ;
plot( b, fb, '*b' ) ;
