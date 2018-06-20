classdef Powell_bsf < FunctionMap
  %
  % Function Powell bad scaled function
  %
  % Detailed reference missing:
  % "Powell,M.J.D.","A hybrid method for nonlinear equations" 1970
  %
  % see also in reference test N.3
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

    function self = Powell_bsf()
      self@FunctionMap(int32(2),int32(2)) ; %
      %exact_solutions        = [];                % no exacts solution provided
      approximated_solutions  = [1.098e-05,9.106]; % Known approximated solution
      guesses                 = [0;1] ;            % one guess
    end

    function F = evalMap(self,x)
      % evaluate the entries (not squared) of the 
      % Powell badly scaled function.
      X1 = x(1) ;
      X2 = x(2) ;
      F = [ (1e04*X1*X2 -1) ; exp(-X1) + exp(-X2) - 1.0001 ] ; % vector of [ f_1(x) ... f_n(x) ] values.
    end

    function J = jacobian( self, x )
      % use analytic jacobian
      self.check_x( x );
      X1 = x(1) ;
      X2 = x(2) ;
      J = [1e04*X2  ,  1e04*X1;   ...
           -exp(-X1),  -exp(-X2)];
    end

    function T = tensor( self, x )
      % use analytic tensor
      self.check_x( x );
      X1 = x(1) ;
      X2 = x(2) ;
      % Create the n-matrices of T

      % D J / D X1
      T1 = [ 0        , 1e04    ;...
             exp(-X1) , 0    ];

      % D J / D X2
      T2 = [ 1e04     , 0      ;...
             0        , exp(-X2) ];
      % Concatenate the n-matrices of T
      % Dimensions = Function map : first D : second D
      T  = cat(3,T1,T2);

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