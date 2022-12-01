addpath('../lib');
addpath('../maps');
addpath('../functions');

%% Trivial test

subplot(2,1,1);
funTest1 = TrivialQuadratic();

fprintf('\n\n == Trivial Quadratic Test == \n\n')
LM1 = MinimizationLevembergMarquardt( funTest1 );
LM1.save_iterate_on();
LM1.setTolerance(1e-6);
LM1.setEpsilon2(1e-10);
LM1.setTau(1e-4);
LM1.setMaxIteration(int32(1000));

x0 = [-1.2;1];
[x_star,converged] = LM1.minimize(x0);

funTest1.contour([-2 2 ],[-2 2 ],50);
hold on;
minCoo = funTest1.exact_solutions;
plot(minCoo(1),minCoo(2),'redx','MarkerSize',10);
LM1.plotIter();

%% Rosenbrock (Banana) test

subplot(2,1,2);
funTest2 = Rosenbrock();

%%

fprintf('\n\n == Rosenbrock Test == \n\n')
LM2 = MinimizationLevembergMarquardt( funTest2 );
LM2.save_iterate_on();
LM2.setTolerance(1e-6);
LM2.setEpsilon2(1e-10);
LM2.setTau(1e-4);
LM2.setMaxIteration(int32(1000));

x0 = [-1.2;1];
[x_star,converged] = LM2.minimize(x0);

funTest2.contour([-2 2 ],[-2 2 ],100);
hold on;
minCoo = funTest2.exact_solutions;
plot(minCoo(1),minCoo(2),'redx','MarkerSize',10);
LM2.plotIter();
