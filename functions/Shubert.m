classdef Shubert < FunctionND
  %
  % The Shubert Function
  %
  %   based on a code of John Burkardt (http://people.sc.fsu.edu/~jburkardt/)
  %
  % Reference:
  %
  %   Richard Brent,
  %   Algorithms for Minimization with Derivatives,
  %   Dover, 2002,
  %   ISBN: 0-486-41998-3,
  %   LC: QA402.5.B74.
  %
  % Author:
  %   Enrico Bertolazzi
  %   Dipartimento di Ingegneria Indutriale
  %   Universita` degli Studi di Trento
  %   email: enrico.bertolazzi@unitn.it
  %

  methods

    function self = Shubert()
      self@FunctionND(int32(2));
      self.exact_solutions = [ 0; 0 ];     % one known solution 
      self.guesses         = [ 0.5; 1.0];
    end

    function f = eval(self,x)
      % evaluate function
      self.check_x(x);
      factor1 = 0.0;
      for i = 1 : 5
        y = i;
        factor1 = factor1 + y * cos ( ( y + 1.0 ) * x(1) + y );
      end
      factor2 = 0.0;
      for i = 1 : 5
        y = i;
        factor2 = factor2 + y * cos ( ( y + 1.0 ) * x(2) + y );
      end
      f = factor1 * factor2;
    end

    function g = grad( self, x )
      % use analitic gradient
      self.check_x(x);
      g = zeros ( 1, 2 );

      factor1 = 0.0;
      df1dx1  = 0.0;
      for i = 1 : 5
        y = i;
        factor1 = factor1 + y * cos ( ( y + 1.0 ) * x(1) + y );
        df1dx1 = df1dx1 - y * ( y + 1.0 ) * sin ( ( y + 1.0 ) * x(1) + y );
      end

      factor2 = 0.0;
      df2dx2 = 0.0;
      for i = 1 : 5
        y = i;
        factor2 = factor2 + y * cos ( ( y + 1.0 ) * x(2) + y );
        df2dx2 = df2dx2 - y * ( y + 1.0 ) * sin ( ( y + 1.0 ) * x(2) + y );
      end

      g(1) = df1dx1 * factor2;
      g(2) = factor1 * df2dx2;
    end

    function h = hessian( self, x )
      % use analitic hessian
      self.check_x(x);
      h = zeros ( 2, 2 );
      factor1 = 0.0;
      df1dx1  = 0.0;
      df1dx11 = 0.0;
      for i = 1 : 5
        y = i;
        factor1 = factor1 + y * cos ( ( y + 1.0 ) * x(1) + y );
        df1dx1 = df1dx1 - y * ( y + 1.0 ) * sin ( ( y + 1.0 ) * x(1) + y );
        df1dx11 = df1dx11 - y * ( y + 1.0 )^2 * cos ( ( y + 1.0 ) * x(1) + y );
      end

      factor2 = 0.0;
      df2dx2 = 0.0;
      df2dx22 = 0.0;
      for i = 1 : 5
        y = i;
        factor2 = factor2 + y * cos ( ( y + 1.0 ) * x(2) + y );
        df2dx2 = df2dx2 - y * ( y + 1.0 ) * sin ( ( y + 1.0 ) * x(2) + y );
        df2dx22 = df2dx22 - y * ( y + 1.0 )^2 * cos ( ( y + 1.0 ) * x(2) + y );
      end

      h(1,1) = df1dx11 * factor2;
      h(1,2) = df1dx1 * df2dx2;
      h(2,1) = df1dx1 * df2dx2;
      h(2,2) = factor1 * df2dx22;
    end
  end
end
