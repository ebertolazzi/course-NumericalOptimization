%
classdef secant_class < handle
  properties  (SetAccess = private, Hidden = true)
    fun
    tol
    max_iter
    verbose
    iter      % iteration done
    flag      % last computation convergence flag
    x_history % saved iteration
  end
  methods
    %
    %  Class contructor
    %  self = this for C++ user
    %
    %  use varargin to initialize object in many ways:
    %
    %  obj = secant_class();
    %  obj = secant_class( tol );
    %  obj = secant_class( tol, max_iter );
    %  obj = secant_class( tol, max_iter, verbose );
    %
    function self = secant_class( varargin )
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
    function setup( self, f )
      self.fun = f;
    end
    %
    % return all the iteration done
    %
    function x_history = get_history( self )
      x_history = self.x_history;
    end
    %
    % Declare the methods solver without inplementation.
    % The implementation is in the file solve.m in the
    % same directory
    %
    x = solve( self, f, x0, x1 )
  end
end
