clear all;
addpath('../lib');
addpath('../functions1D');

% define the object function fs, fc
fs = SinFun();

x = 0:2*pi/1000:2*pi;

subplot( 3, 1, 1);
hold off;
disp('SinFun');
plot( x, fs.eval(x), '-k', 'Linewidth', 2 );
title('SinFun');
hold on;

subplot( 3, 1, 2);
hold off;
disp('SinFun derivative with finite difference');
plot( x, fs.eval_D(x), '-r', 'Linewidth', 3 );
hold on;
plot( x, fs.FD_eval_D(x), '--b', 'Linewidth', 2 );
title('SinFun''');

subplot( 3, 1, 3);
hold off;
disp('SinFun derivative with finite difference');
plot( x, fs.eval_DD(x), '-r', 'Linewidth', 3 );
hold on;
plot( x, fs.FD_eval_DD(x), '--b', 'Linewidth', 2 );
title('SinFun''''');
hold on;
