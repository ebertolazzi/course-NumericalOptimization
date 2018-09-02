classdef Barrier1 < FunctionND
  % Author:
  %   Enrico Bertolazzi
  %   Dipartimento di Ingegneria Indutriale
  %   Universita` degli Studi di Trento
  %   email: enrico.bertolazzi@unitn.it
  %
  properties (SetAccess = private, Hidden = true)
    epsi;
  end

  methods

    function self = Barrier1()
      self@FunctionND(int32(2));
      self.exact_solutions = []; % no known solution
      self.guesses         = [ 0.9999; 0];
      self.epsi            = 1e-10;
    end

    function f = eval(self,xx)
      % evaluate function
      self.check_x(xx);
      x = xx(1);
      y = xx(2);
      c = 1-x^2-y^2;
      if c <= 0
        f = inf;
      else
        f = x+2*y+self.epsi/c;
      end
    end

    function g = grad( self, xx )
      % use analitic gradient
      self.check_x(xx);
      x = xx(1);
      y = xx(2);
      c = 1-x^2-y^2;
      g = [1,2];
      if c > 0
        g = g + (2*self.epsi/c^2)*[x,y];
      end
    end

    function H = hessian( self, xx )
      % use analitic hessian
      self.check_x(xx);
      x = xx(1);
      y = xx(2);
      c = 1-x^2-y^2;
      H = zeros( 2, 2 );
      if c > 0
        H = (2*self.epsi/c^3)*[ y^2-3*x^2-1, -4*x*y;...
                           -4*x*y,       x^2-3*y^2-1 ];
      end
    end
  end
end
