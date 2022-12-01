%
classdef halley_class < handle
  properties  (SetAccess = private, Hidden = true)
    fun
    Dfun
    DDfun
    tol
    max_iter
    verbose
    iter      % iteration done
    flag      % last computation convergence flag
    x_last    % last computed value
  end
  methods
    %
    %  Class contructor
    %  self = this for C++ user
    %
    %  use varargin to initialize object in many ways:
    %
    %  obj = halley_class();
    %  obj = halley_class( tol );
    %  obj = halley_class( tol, max_iter );
    %  obj = halley_class( tol, max_iter, verbose );
    %
    function self = halley_class( varargin )
      % set default values
      self.tol      = 1e-8;
      self.max_iter = 100;
      self.verbose  = 'iter';
      if nargin > 0
        %
        % in a true class require to check the input
        % to be correct, 1 real number > 0 etc.
        %
        self.tol = varargin{1};
      end
      if nargin > 1
        self.max_iter = varargin{2};
      end
      if nargin > 2
        self.verbose = varargin{3};
      end
    end
    %
    function set_tolerance( self, new_tol )
      self.tol = new_tol;
    end
    %
    function set_max_iter( self, max_iter )
      self.max_iter = max_iter;
    end
    %
    function set_verbose( self, verbose )
      self.verbose = verbose;
    end
    %
    function iter = get_iter( self )
      iter = self.iter;
    end
    %
    % store the function f and its derivative Df
    % in the class
    %
    function setup( self, f, Df, DDf )
      self.fun   = f;
      self.Dfun  = Df;
      self.DDfun = DDf;
    end
    %
    % Declare the methods solver without inplementation.
    % The implementation is in the file solve.m in the
    % same directory
    %
    x = solve( self, f, Df, DDf, x0 )
  end
end
