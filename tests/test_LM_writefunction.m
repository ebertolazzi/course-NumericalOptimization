
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
               ' you the minimum of:\nExample: ExtendedRosenbrock( 10 )\n\n'],'s');

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

disp('Minimum found, it is:      x = ');
disp(x_star);

% display the solutions available from the class

if ~isempty(funTest1.exact_solutions)
    disp('Exact solutions:');
    disp(funTest1.exact_solutions);
elseif ~isempty(funTest1.approximated_solutions)
    disp('Approximated solutions:');
    disp(funTest1.approximated_solutions);
else
    disp('It was not possible to find given solutions in the function class');
end








