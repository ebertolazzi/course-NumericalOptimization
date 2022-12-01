clear all;
addpath('../lib');
addpath('../functions1D');

FUN = testFun2();
%FUN = testFun3();
%FUN = testFun4();

% initialize Golden Search object m
m = LinesearchGoldenSection();
% attach function fs to the miminization object
m.setFunction(FUN);

aa = 0;
bb = 10;
x  = aa:(bb-aa)/1000:bb;

subplot( 2, 1, 1);
hold off;
plot( x, FUN.eval(x), '-k', 'Linewidth', 1 );
hold on;

disp('do golden search');
%%
[ LO, HI, ierr ] = m.ForwardBackward( bb/3 );
fprintf('ierr = %d\n',ierr);
disp(struct2table([LO,HI]));

plot( x, m.ArmijoSlope(x), '-g', 'Linewidth', 1 );
title('testFun');

x = LO.alpha:(HI.alpha-LO.alpha)/1000:HI.alpha;
plot( x, FUN.eval(x), '-r', 'Linewidth', 3 );

[a,b] = m.minimize(LO.alpha, HI.alpha);
%%
subplot( 2, 1, 2);
hold off;
x = aa:(bb-aa)/1000:bb;
plot( x, FUN.eval(x), '-r', 'Linewidth', 1 );
hold on;
plot( x, m.ArmijoSlope(x), '-g', 'Linewidth', 1 );

x = a:(b-a)/1000:b;
plot( x, FUN.eval(x), '-b', 'Linewidth', 2 );
hold on;
plot( (a+b)/2, FUN.eval((a+b)/2), 'ok' );

title('Golden Search');
