
% TEST (only 2D)

addpath(genpath('/Users/giammarco/GoogleDriveUNITN/PhDall/corsi/BertoProject/optimization/lib'));
addpath(genpath('/Users/giammarco/GoogleDriveUNITN/PhDall/corsi/BertoProject/optimization/functions'));



%% Trivial test

funTest1 = TrivialQuadratic();

%%

fprintf(1,'\n\n == Trivial Quadratic Test == \n\n')

LM1    = MinimizationLM( funTest1 );

LM1.debug_on;

LM1.setTolerance(1e-6);
LM1.setEpsilon2(1e-10);

LM1.setTau(1e-4);
LM1.setMaxIteration(int32(1000));

x0 = [-1.2;1];

[x_star,converged] = LM1.minimize(x0);

% draw the contoLM1ur of the function (only 2D)

xlim = [-2 2 ];
ylim = [-2 2 ];

[X,Y] = meshgrid(xlim(1):0.01:xlim(2),ylim(1):0.01:ylim(2));

Z = zeros(size(X));

for i = 1:size(X,2)
    for j = 1:size(Y,1)
        
        Z(j,i) = funTest1.eval([ X(1,i) ; Y(j,1) ]);
        
    end
end

%%

h = figure(1);
title('Trivial quadratic function')
contour(X,Y,Z,logspace(log10(0.001),log10(max(max(Z))),100));
axis([xlim ylim]);
hold on;
minCoo = funTest1.exact_solutions;
plot(minCoo(1),minCoo(2),'redx','MarkerSize',10);
hold on;
x = LM1.x_history(1,:);
y = LM1.x_history(2,:);
plot(x,y,'blue','Marker','o','LineWidth',1,'color',[1 0.7 0.1]);
hold off;
clear x y;

movegui(h,'northwest')

testPassedDisplay( x_star , funTest1 , 10^-4 );






%% Resenbrock (Banana) test

funTest2 = Rosenbrock();

%%

fprintf(1,'\n\n == Rosenbrock Test == \n\n')

LM2    = MinimizationLM( funTest2 );

LM2.debug_on;

LM2.setTolerance(1e-6);
LM2.setEpsilon2(1e-10);

LM2.setTau(1e-4);
LM2.setMaxIteration(int32(1000));

x0 = [-1.2;1];

[x_star,converged] = LM2.minimize(x0);

% draw the contofunTest2ur of the function (only 2D)

xlim = [-2 2 ];
ylim = [-2 2 ];

[X,Y] = meshgrid(xlim(1):0.01:xlim(2),ylim(1):0.01:ylim(2));

Z = zeros(size(X));

for i = 1:size(X,2)
    for j = 1:size(Y,1)
        
        Z(j,i) = funTest2.eval([ X(1,i) ; Y(j,1) ]);
        
    end
end

%%

h = figure(2);
title('Banana function (Rosenbrock)')
contour(X,Y,Z,logspace(log10(0.001),log10(max(max(Z))),100));
axis([xlim ylim]);
hold on;
minCoo = funTest2.exact_solutions;
plot(minCoo(1),minCoo(2),'redx','MarkerSize',10);
hold on;
x = LM2.x_history(1,:);
y = LM2.x_history(2,:);
plot(x,y,'blue','Marker','o','LineWidth',1,'color',[1 0.7 0.1]);
hold off;
clear x y;

movegui(h,'northeast')



testPassedDisplay( x_star , funTest2 , 10^-4 );


function testPassedDisplay( x_star , funTest , passed_test_sol )
	if ( norm( x_star - funTest.exact_solutions() ) < abs( passed_test_sol ) )
		fprintf(1,'test passed and the minimum is in [ %3.3g %3.3g ] \n' , x_star(1) , x_star(2) );
	else
		fprintf(1,'test not passed, fake minimum is in [ %3.3g %3.3g ], come back the next session \n' , x_star(1) , x_star(2) );
	end
end

function testPassedApproxDisplay( x_star , funTest , passed_test_sol )
	if ( norm( x_star - funTest.exact_solutions() ) < abs( passed_test_sol ) )
		fprintf(1,'test passed and the minimum is in [ %3.3g %3.3g ] \n' , x_star(1) , x_star(2) );
	else
		fprintf(1,'test not passed, fake minimum is in [ %3.3g %3.3g ], come back the next session \n' , x_star(1) , x_star(2) );
	end
end

