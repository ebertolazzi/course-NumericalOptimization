addpath('../functions');
addpath('../lib');

% get the linesearch object and initialize it to the cut function FUN1D
LFB = LinesearchMoreThuente();
%LFB = LinesearchGoldenSection();

% cut the function near the barrier and plot
%FUN1D = ls1();
%FUN1D = ls2();
%FUN1D = ls3();
%FUN1D = ls4();
%FUN1D = ls5();%%%%%%%%%%
%FUN1D = ls6();
%FUN1D = ls7();
%FUN1D = ls8();
FUN1D = ls9();
LFB.setFunction( FUN1D );

range = FUN1D.getRange();
t     = (0:range(2)/1000:range(2));
y     = zeros(size(t));
y1    = zeros(size(t));
z     = zeros(size(t));

f0  = FUN1D.eval(0)
Df0 = FUN1D.eval_D(0)
z   = f0 + LFB.c1*t*Df0;

for k=1:length(t)
  y(k)  = FUN1D.eval(t(k));
  y1(k) = FUN1D.eval(-t(k)/5);
end
hold off
plot( t, y, '-r', 'Linewidth', 3 );
hold on
%plot( -t/5, y1, '-b', 'Linewidth', 1 );
plot( t, z, '-g', 'Linewidth', 1 );

for alpha_guess=[1e-6,0.1,1,10,1000]
%for alpha_guess=[1]
  alpha = LFB.search( alpha_guess );
  plot( alpha, FUN1D.eval(alpha), 'ob', 'Linewidth', 3 );
  alpha
end

