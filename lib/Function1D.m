classdef (Abstract) Function1D < handle
  % This class implement two methods (apart from the constructor), which
  % evaluate the first and second derivative of a 1-D function. The
  % derivation is performed with the finite difference rule, so it's not
  % analytical but is an approximation.
  % The method here proposed is not the classical one direction incremental
  % ratio, but is the bi-directional one (more accurate)

  methods (Abstract)
    %
    % Define the abstract functions.
    % eval, eval_D, eval_DD, eval_FG and eval_FGH
    % will be defined in the derived class
    y          = eval( self, x )
    Dy         = eval_D( self, x )
    DDy        = eval_DD( self, x )
    [y,Dy]     = eval_FG( self, x )
    [y,Dy,DDy] = eval_FGH( self, x )
  end

  methods
    % Finite difference approximation
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function Dy = FD_eval_D( self, x )
      % Finite difference approximation of the first derivative
      h  = max(1,abs(x))*sqrt(eps); % "eps" is the "machine epsilon precision"
      fp = self.eval(x+h);
      fm = self.eval(x-h);
      Dy = (fp-fm)./(2*h);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function DDy = FD_eval_DD( self, x )
      % Finite difference approximation of the second derivative
      h   = max(1,abs(x))*(eps^(1/3));
      fp  = self.eval(x+h);
      fm  = self.eval(x-h);
      fc  = self.eval(x);
      DDy = (fp+fm-2*fc)./(h.^2);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
