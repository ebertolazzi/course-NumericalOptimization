clear all;
addpath('../functions1D');
addpath('../lib');

LS1 = ls1();
LS2 = ls2();
LS3 = ls3();
LS4 = ls4();
LS5 = ls5();
LS6 = ls6();
LS7 = ls7();
LS8 = ls8();
LS9 = ls9();

subplot(3,3,1);
t = 0:0.001:6;
plot( t, LS1.eval(t), 'Linewidth', 3 );

subplot(3,3,2);
t = 0:0.001:2;
plot( t, LS2.eval(t), 'Linewidth', 3 );

subplot(3,3,3);
t = 0:0.001:3;
plot( t, LS3.eval(t), 'Linewidth', 3 );

subplot(3,3,4);
t = 0:0.001:1;
plot( t, LS4.eval(t), 'Linewidth', 3 );

subplot(3,3,5);
t = 0:0.001:1;
plot( t, LS5.eval(t), 'Linewidth', 3 );

subplot(3,3,6);
t = 0:0.001:1;
plot( t, LS6.eval(t), 'Linewidth', 3 );


subplot(3,3,7);
t = 0:0.001:2;
plot( t, LS7.eval(t), 'Linewidth', 3 );

subplot(3,3,8);
t = 0:0.001:0.4;
plot( t, LS8.eval(t), 'Linewidth', 3 );

subplot(3,3,9);
t = 0:0.001:10;
plot( t, LS9.eval(t), 'Linewidth', 3 );
