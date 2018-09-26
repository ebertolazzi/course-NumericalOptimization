classdef GaussianFunction < FunctionMap
  %
  % Gaussian Function N = 3 , M = 15;
  %
  % unpublished
  % 
  % see also number 9 in:
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

  properties( SetAccess = private ,  Hidden = true )
    yi_i % table for y_i given i ....... see main reference: More, Jorge J. and Garbow, Burton S. and Hillstrom, Kenneth E.,
         % Testing Unconstrained Optimization Software. (function number 9)
    yi_half
  end

  methods

    function self = GaussianFunction()
      % GaussianFunction()...............N = 3 ; M = 15;
      self@FunctionMap(int32(3),int32(15)) ;        % call superclass constructor (initialize M)
      %exact_solutions        = [];                % no exacts solution provided
      self.approximated_solutions =  1.12793 * 10^-8   ;                 % approximated solution provided only if M == 10
      self.guesses                = [ 0.4 , 1 , 0     ];                 % one guess

      self.yi_half = [ 0.0009 0.0044 0.0175 0.0540 0.1295 0.2420 0.3521 0.3989 ];
      self.yi_i    = [ self.yi_half fliplr(self.yi_half(1:end-1)) ].'; % it's a symmetric vector w.r.t to i = 8;
                                                             % it's a column vector for consistency with Jacobian
                                                             % actually it is a gaussian -.-      
    end

    function F = evalMap(self,x)
      % evaluate the entries (not squared) of the function.
      X1 = x(1) ;
      X2 = x(2) ;
      X3 = x(3) ;
      i    = (1:self.M).'  ; % create indexes for yi_i
      ti_i = ( 8 - i ) ./ 2; % ti_i for the function
      F    = ( X1 .* exp( ( - X2 .* ( ti_i - X3 ) .^ 2 ) ./ 2) ) - self.yi_i ; % vector of [ f_1(x) ... f_n(x) ] values.
                                                                            % i is automatically spanned
    end

    function J = jacobian( self, x )
      % use analytic jacobian
      self.check_x( x );
      X1 = x(1) ;
      X2 = x(2) ;
      X3 = x(3) ;
      i    = (1:self.M).'  ; % create indexes for yi_i
      ti_i = ( 8 - i ) ./ 2; % ti_i for the function
      J  = [ exp(-X2 .* (ti_i - X3) .^ 2 / 2),...
                 -X1 .* (ti_i - X3) .^ 2  .* exp( -X2 .* (ti_i - X3) .^ 2 ./ 2) ./ 2, ...
                  X1 .* X2 .* (ti_i - X3) .* exp( -X2 .* (ti_i - X3) .^ 2 ./ 2)];
    end

    function T = tensor( self, x )
      % use analytic tensor
      self.check_x( x );
      X1 = x(1) ;
      X2 = x(2) ;
      X3 = x(3) ;
      % Create the n-matrices of T

      % D J / D X1
      T1 = [zeros(self.M,1), - ( ti_i - X3) .^ 2 .* exp(-X2 .* (ti_i - X3) .^ 2 / 2) / 2, X2 .* (ti_i - X3) .* exp(-X2 .* (ti_i - X3) .^ 2 ./ 2)];

      % D J / D X2
      T2 = [-(ti_i - X3) .^ 2 .* exp(-X2 .* (ti_i - X3) .^ 2 / 2) / 2,...
       X1 .* (ti_i - X3) .^ 4 .* exp(-X2 .* (ti_i - X3) .^ 2 / 2) / 4,...
        X1 .* (ti_i - X3) .* exp(-X2 .* (ti_i - X3) .^ 2 ./ 0.2e1) - X1 .* X2 .* (ti_i - X3) .^ 3 .* exp(-X2 .* (ti_i - X3) .^ 2 ./ 2) ./ 2];

      % D J / D X3
      T3 = [ X2 .* (ti_i - X3) .* exp(-X2 .* (ti_i - X3) .^ 2 ./ 2),...
       X1 .* (ti_i - X3) .* exp(-X2 .* (ti_i - X3) .^ 2 / 2) - X1 .* X2 .* (ti_i - X3) .^ 3 .* exp(-X2 .* (ti_i - X3) .^ 2 ./ 2) ./ 2,...
        -X1 .* X2 .* exp(-X2 .* (ti_i - X3) .^ 2 ./ 2) + X1 .* X2 .^ 2 .* (ti_i - X3) .^ 2 .* exp(-X2 .* (ti_i - X3) .^ 2 ./ 2)];



      % Concatenate the n-matrices of T
      % Dimensions = MxNxN
      T  = cat(3,T1,T2,T3);

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