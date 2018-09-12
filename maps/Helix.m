classdef Helix < FunctionMap
  %
  % The Helix function.
  %
  % Reference:
  %   Fletcher, R. and Powell, M. J. D.
  %   A Rapidly Convergent Descent Method for Minimization
  %   The Computer Journal, vol 6, n.2, pp. 163-168, 1963,
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
  methods (Hidden = true)
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = f_theta(self,X,Y)
      res = atan2(Y,X)/(2*pi);
      if X < 0
        res = res + 0.5;
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = f_theta_1(self,X,Y)
      res = -Y/(X^2+Y^2);
      res = res/(2*pi);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = f_theta_2(self,X,Y)
      res = X/(X^2+Y^2);
      res = res/(2*pi);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = f_theta_1_1(self,X,Y)
      res = 2*X*Y/(X^2+Y^2)^2;
      res = res/(2*pi);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = f_theta_1_2(self,X,Y)
      res = (Y^2-X^2)/(X^2+Y^2)^2;
      res = res/(2*pi);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = f_theta_2_2(self,X,Y)
      res = -2*X*Y/(X^2+Y^2)^2;
      res = res/(2*pi);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = hypot_1(self,x,y)
      res = x/hypot(x,y);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = hypot_2(self,x,y)
      res = y/hypot(x,y);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = hypot_1_1(self,x,y)
      res = y^2/hypot(x,y)^(3/2);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = hypot_1_2(self,x,y)
      res = -x*y/hypot(x,y)^(3/2);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = hypot_2_2(self,x,y)
      res = x^2/hypot(x,y)^(3/2);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end

  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Helix()
      self@FunctionMap(int32(3),int32(3));
      self.exact_solutions = [ 1; 0; 0 ];
      self.guesses         = [ -1; 0; 0 ];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function F = evalMap(self,x)
      self.check_x(x);
      X = x(1);
      Y = x(2);
      Z = x(3);
      F = [ 10*( Z - 10*self.f_theta(X,Y) ); 10*( hypot(X,Y)  - 1 ); Z ];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function J = jacobian( self, x )
      % use analitic jacobian
      self.check_x(x);
      X = x(1);
      Y = x(2);
      Z = x(3);
      J = [ -100*self.f_theta_1(X,Y), -100*self.f_theta_2(X,Y), 10; ...
            10*self.hypot_1(X,Y), 10*self.hypot_2(X,Y), 0; ...
            0, 0, 1 ];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function T = tensor( self, x )
      % use analitic tensor of second derivative
      X = x(1);
      Y = x(2);
      Z = x(3);
      T        = zeros(3,3,3);
      T(1,:,:) = [ -100*self.f_theta_1_1(X,Y), -100*self.f_theta_1_2(X,Y), 0; ...
                   -100*self.f_theta_1_2(X,Y), -100*self.f_theta_2_2(X,Y), 0; ...
                   0, 0, 0 ];
      tmp = 10/(X^2+Y^2)^(3/2);
      T(2,:,:) = tmp*[ Y^2, -X*Y, 0; -X*Y, X*X, 0; 0, 0, 0 ];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
