classdef Powell_bsf < FunctionMap
  %
  % Function Powell bad scaled function
  %
  % ORIGINAL REFERENCE MISSING
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
  % Author: Giammarco Valenti -UNCOMPLETE VERSION- 0.0
  %
  methods

    function self = Powell_bsf()
      self@FunctionMap(int32(2),int32(2)) ;
      %exact_solutions        = [];                % no exacts solution provided
      approximated_solutions  = [1.098e-05,9.106]; % Known approximated solution
      guesses                 = [0;1] ;            % one guess
    end

    function F = evalMap(self,x)
      % evaluate the entries (not squared) of the 
      % Powell badly scaled function.
      X = x(1) ;
      Y = x(2) ;
      F = [ (1e04*X*Y -1) ; exp(-X) + exp(-Y) - 1.0001 ] ; % vector of [ f_1(x) ... f_n(x) ] values.
    end

    function J = jacobian( self, x )
      % use numeric jacobian
      self.check_x( x );
      J = self.FD_jacobian( x ); 
    end

    function T = tensor( self, x )
      % use numeric tensor of second derivative
      self.check_x( x );
      T = self.FD_tensor( x );
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

% Coded by Giammarco Valenti - 2018-01-07 - version 0.0 - ☠☠☠ Beta
% /// KEEP IT SIMPLE! \\\