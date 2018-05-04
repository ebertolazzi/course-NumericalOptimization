classdef JennrichAndSampson < FunctionMap
  %
  % Jennrich And Sampson function , M as argument
  %
  % JENNRICH, R.I., AND SAMPSON,P.F. Applicatmn of stepwise
  % regression to nonlinear estimatmn. Technometrtcs 10 (1968), 63-72.
  %
  % see also in reference test N.6
  %
  % @article{More:1981,
  %   author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and Hillstrom, Kenneth E.},
  %   title   = {Testing Unconstrained Optimization Software},
  %   journal = {ACM Trans. Math. Softw.},
  %   volume  = {7},
  %   number  = {1},
  %   year    = {1981},
  %   pages   = {17--41},
  %   doi     = {10.1145/355934.355936}
  % }
  %
  % Author: Giammarco Valenti - University of Trento
  % -> final debug required
  methods

    function self = JennrichAndSampson( M )
      % JennrichAndSampson( M )...............N = 2 ; M = M;
      self@FunctionMap(int32(2),int32(M)) ;        % call superclass constructor (initialize M)
      %exact_solutions        = [];                % no exacts solution provided
      if M == 10
        approximated_solutions = 0.2578.*[1 1];         % approximated solution provided only if M == 10
      end
      guesses                 = [0.3 0.4].' ;            % one guess
    end

    function F = evalMap(self,x)
      % evaluate the entries (not squared) of the 
      % Powell badly scaled function.
      X1 = x(1) ;
      X2 = x(2) ;
      i  = (1:self.M).'; % column vector required
      F  = 2 + 2.*i - ( exp( i.*X1 ) + exp( i.*X2 ) ); % vector of [ f_1(x) ... f_n(x) ] values.
    end

    function J = jacobian( self, x )
      % use analytic jacobian
      self.check_x( x );
      X1 = x(1) ;
      X2 = x(2) ;
      i  = (1:self.M).'; % column vector required
      J  = [ -i .* exp(i .* X1) , -i .* exp(i .* X2) ];
    end

    function T = tensor( self, x )
      % use analytic tensor
      self.check_x( x );
      X1 = x(1) ;
      X2 = x(2) ;
      i  = (1:self.M).'; % column vector required
      % Create the n-matrices of T

      % D J / D X1
      T1 = [ -i .^ 2 .* exp( i .* X1 ) , zeros(self.M,1)        ];

      % D J / D X2
      T2 = [ zeros(self.M,1)           , -i .^ 2 .* exp( i .* X2 ) ];

      % Concatenate the n-matrices of T
      % Dimensions = MxNxN
      T  = cat(3,T1,T2);

    end

    % For tensor and jacobian a maple file is available: ask to the author if needed

  end
end

%  #######################################################
%  #  _______   __  _____ ______ _______ _ ____   _____  #
%  #         \ |  |   ___|   ___|__   __| |    \ |       #
%  #          \|  |  __| |  __|    | |  | |     \|       #
%  #       |\     | |____| |       | |  | |   |\         #
%  #  _____| \____| _____|_|       |_|  |_|___| \______  #
%  #                                                     #
%  #######################################################