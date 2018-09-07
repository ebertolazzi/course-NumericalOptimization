addpath('../functions');
addpath('../lib');

% get a 2d barrier function
FUN2D = Barrier1();
%x0    = [0.99999999999; 0];
x0    = [0.999993944545; 0];
d     = [ -1; 0*(-1)];

x0 = [-0.447147077055711;-0.894293265304308];
g  = FUN2D.grad( x0 );
d  = -g.'/norm(g);

% cut the function near the barrier and plot
FUN1D = Function1Dcut( FUN2D, x0, d );
FUN1D.no_FD_D();
tmax = 1e-10;

subplot(2,1,1);
t  = (0:tmax/1000:tmax);
y  = zeros(size(t));
y1 = zeros(size(t));
for k=1:length(t)
  y(k)  = FUN1D.eval(t(k))-FUN1D.eval(0);
  y1(k) = FUN1D.eval(-t(k))-FUN1D.eval(0);
end
hold off
plot( t, y, '-r', 'Linewidth', 3 );
hold on
plot( -t, y1, '-b', 'Linewidth', 1 );

% get the linesearch object and initialize it to the cut function FUN1D
%LFB = LinesearchForwardBackward('check LS');
LFB = LinesearchMoreThuente();
LFB.setFunction( FUN1D );

alpha_guess = 2.1;
alpha       = LFB.search( alpha_guess );

alpha
LFB.printInfo();

subplot(2,1,2);
t   = (0:alpha/1000:alpha);
y   = zeros(size(t));
z   = zeros(size(t));
f0  = FUN1D.eval(0)
Df0 = 0.2*FUN1D.eval_D(0)
for k=1:length(t)
  y(k) = min([2,FUN1D.eval(t(k))]);
  z(k) = f0+t(k)*Df0;
end
hold off
plot( t, y, '-r', 'Linewidth', 1 );
hold on
plot( t, z, '-b', 'Linewidth', 3 );
hold off
