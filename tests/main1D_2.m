clear all;
addpath('../lib');

c1  = 0.1 ;
c2  = 0.5 ;
f1D = testFun() ;

LS = LinesearchArmijo();
LS.setC1C2( c1, c2 );
LS.setFunction( f1D );

x = 0:pi/1000:pi ;

%[alpha,ok] = LS.search( 0.9*pi ) ;
[alpha0,alpha1,ok] = LS.ForwardBackward( 1e-10*pi );


fprintf('alpha = %g, ok = %g\n',alpha0,ok);

hold off ;
plot( x, f1D.eval(x) ) ;
hold on ;
f0  = f1D.eval(0) ;
Df0 = f1D.eval_D(0) ;
r   = @(xx) f0+c1*xx*Df0 ;
plot( x, r(x), '-r' ) ;
plot( alpha0, r(alpha0), 'ob' ) ;
plot( alpha0, f1D.eval(alpha0), 'ok' ) ;
plot( alpha1, r(alpha1), 'xb' ) ;
plot( alpha1, f1D.eval(alpha1), 'xk' ) ;
