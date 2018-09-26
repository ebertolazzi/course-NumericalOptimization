classdef ExtendedRosenbrock < FunctionMap
  %
  %
  % @article{spedicato1975computational,
  %   title={Computational experience with quasi-Newton algorithms for minimization problems of moderately large size},
  %   author={Spedicato, E},
  %   year={1975}
  % }
  %
  % see also number 21 in:
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
  methods

    function self = ExtendedRosenbrock( M )
      % Extendedrosenbrock( M ) M=N ... has to be even!
      if ( mod(M,2) ~= 0 | M <= 0 )
        error('ExtendedRosenbrock:: Bad constructor usage: M ( M = N ) has to be EVEN and greater than 0')
      end
      N = M; % The same
      self@FunctionMap( int32(N) , int32(M) );          % call superclass constructor M = N
      self.exact_solutions        = ones( self.N , 1 );      % no exacts solution provided
      %self.approximated_solutions = [];                % approximated solution provided only 
      self.guesses                = ones( self.N , 1 );
      self.guesses(1:2:self.N)         = -1.2         ;      % one guess
    end

    function F = evalMap( self , x )
      % evaluate the entries (not squared) of the function (row vector)
      odd   = 1:2:self.N; % M = self.N
      even  = 2:2:self.N;

      F         = zeros( self.M , 1 ); 
      F( odd )  = 10.*( x( even ) - x( odd ).^2 );
      F( even ) = 1 - x( odd );

      % F is the vector of [ f_1(x) ... f_n(x) ] values.
    end

    function J = jacobian( self, x )
      % use analytic jacobian
      self.check_x( x );
      odd   = 1:2:self.N; % M = self.N
      even  = 2:2:self.N;

      J  = zeros(self.N);

      N = self.N;

      % MATLAB CHALLENGE! square-by-square diagonal matrix

      J( 1   : 2*N+2 : N*N ) = -20 * x( odd );
      J( 2   : 2*N+2 : N*N ) = -1;
      J( N+1 : 2*N+2 : N*N ) = 10;

      J = J.';


    end

    function T = tensor( self, x )
      % use analytic tensor
      self.check_x( x );

      N = self.N;

      T  = zeros(N,N,N);

      % MATLAB CHALLENGE! uno-ogni-tanto non zero pseudo-diagonal tensor

      T(1 : (2*N*N+2*N+2) : N*N*N ) = -20;



    end
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