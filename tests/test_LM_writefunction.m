
% TEST LM choose the function

%  #######################################################
%  #  _______   __  _____ ______ _______ _ ____   _____  #
%  #         \ |  |   ___|   ___|__   __| |    \ |       #
%  #          \|  |  __| |  __|    | |  | |     \|       #
%  #       |\     | |____| |       | |  | |   |\         #
%  #  _____| \____| _____|_|       |_|  |_|___| \______  #
%  #                                                     #
%  #######################################################

% Giammarco Valenti

addpath(genpath('../lib'));
addpath(genpath('../functions'));

commando = input(['\nWrite the "call" to the function you want Levendberg Marquart to provide' ...
               'you the minimum of:\nExample: BoxThreeDimensionalFunction( 10 )\n\n'],'s');

eval([ 'funTest1 =' commando ]);


fprintf(1,[ '\n\n == ' commando ' == \n\n' ])

LM1    = MinimizationLM( funTest1 );

LM1.debug_on;

LM1.setTolerance(1e-6);
LM1.setEpsilon2(1e-10);

LM1.setTau(1e-2);
LM1.setMaxIteration(int32(100));

x0 = funTest1.guesses;

[x_star,converged] = LM1.minimize(x0); % Minimization

x_ex = funTest1.exact_solutions; % Extract exact solutions

disp('Minimum found, it is:      x = ');
disp(x_star);

try 
    disp('Exact solutions were:      x = '),
    disp(x_ex);
catch
    disp('It was not possible to find exact solutions in the class function');
end








