classdef (Abstract) Minimization1D < handle

  properties (SetAccess = private, Hidden = true)
    funObj   % classe tipo fiunzione 1D
    tol      % tolleranza per |b-a|
    max_iter % massimo numero iterate ammesse
  end

  methods
    function self = Minimization1D( varargin )
      % funObj will be initialized using setFunction
      % varargin{1}  tol
      % varargin{2}  max_iter
      self.tol      = 1e-3 ;
      self.max_iter = 10 ;
      self.funObj   = @(x) error('you must call setFunction before use Minimization1D') ;
      if nargin > 0 ; self.tol      = varargin{1} ; end
      if nargin > 1 ; self.max_iter = varargin{2} ; end
    end

    function setFunction( self, f )
      self.funObj = f ;
    end
  end

  %methods (Abstract)
  %  y = minimize( self, a, b )
  %end

end
