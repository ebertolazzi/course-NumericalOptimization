classdef MinimizationGradientMethod < handle

  properties (SetAccess = private, Hidden = true)
    funND      % class function with expected methods
               % funND.eval(x)
               % funND.grad(x)
               % funND.hessian(x) (UNUSED)
    linesearch % quasi minimizzation with expected method
               % linesearch.setFunction(fun)
               % [alpha,ok] = linesearch.search(alpha_guess)
    tol        % tolleranza per |grad|
    max_iter   % massimo numero iterate ammesse
    debug_on   % if true do some additional printing
    FD_D 
  end

  methods
    function self = MinimizationGradientMethod( linesearch, varargin )
      % varargin{1} tol
      % varargin{2} max_iter
      self.tol        = 1e-3 ;
      self.max_iter   = 100 ;
      self.linesearch = linesearch ;
      self.debug_on   = false ;
      self.FD_D       = true ;
      
      self.funND      = @(x) error('you must call setFunction before use Minimization1D') ;
      if nargin > 1 ; self.tol      = varargin{1} ; end
      if nargin > 2 ; self.max_iter = varargin{2} ; end
      if nargin > 3 ; self.debug_on = varargin{3} ; end
    end

    function setFunction( self, f )
      self.funND = f ;
    end    

    function self = setDebug( self )
      self.debug_on = true ;
    end

    function self = setNoDebug( self )
      self.debug_on = false ;
    end

    function self = use_FD_D( self )
      self.FD_D = true ;
    end

    function self = no_FD_D( self )
      self.FD_D = false ;
    end
    
    function [x1,alpha] = step1D( self, x0, d, alpha_guess )
      %
      % build the 1D function along the search direction
      fcut = Function1Dcut( self.funND, x0, d );
      %
      % set analitic gradient if necessary
      if self.FD_D
        fcut.use_FD_D() ;
      else
        fcut.no_FD_D() ;
      end
      %
      % do a 1D minimization
      self.linesearch.setFunction( fcut ) ;
      %
      % search an interval for minimization
      [alpha,ok] = self.linesearch.search( alpha_guess ) ;
      %
      % check error
      if ~ok
        error('MinimizationGradientMethod:step1D, linesearch failed\n');
      end
      %
      % advance
      x1 = x0 + alpha * d ;
    end
    
    function [xs,converged] = minimize( self, x0 )
      xs    = x0 ;
      alpha = 1 ;
      converged = false ;
      for iter=1:self.max_iter
        %
        % find search direction = - gradient of the function
        d  = -self.funND.grad( xs ).' ;
        %
        % check if converged
        converged = norm(d,inf) < self.tol ;
        if converged ; break ; end 
        %
        % only for debug
        if self.debug_on
          fprintf(1,'iter = %5d ||grad f|| = %14g ...', iter, norm(d,Inf) ) ;
        end
        %
        % minimize along search direction
        [xs,alpha] = self.step1D( xs, d, alpha ) ;
        %
        % only for debug
        if self.debug_on
          fprintf(1,' alpha = %g\n', alpha ) ;
        end
        %
      end
    end

  end
end
