addpath('../functions');
addpath('../lib');
addpath('../maps');

funTest1 = Wood();

TR1 = MinimizationTrustRegion( funTest1 );
TR1.save_iterate_on();
TR1.setTolerance(1e-6);
%TR1.setEpsilon2(1e-10); metti poi tutti i tuoi
%TR1.setTau(1e-10); metti poi anche la scelta dell'algoritmo
TR1.setMaxIteration(int32(100));
TR1.selectByNumber(3);

x0 = funTest1.guesses;

[x_star,converged] = TR1.minimize(x0); % Minimization

x_ex  = funTest1.exact_solutions; % Extract exact solutions

if converged
    fprintf('\nMinimum found');
   % fprintf('\nMinimum found, it is:      x = [ %2.2g %2.2g ]''\n', x_star(1),x_star(2) );
   % fprintf('\nExact solutions were:      x = [ %2.2g %2.2g ]'''  , x_ex(1),x_ex(2));
else
    fprintf('\nNot converged');
end    