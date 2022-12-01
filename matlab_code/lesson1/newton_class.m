%
classdef newton_class < handle
  properties
    fun
    Dfun
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
    %  obj = newton_class();
    %  obj = newton_class( tol );
    %  obj = newton_class( tol, max_iter );
    %  obj = newton_class( tol, max_iter, verbose );
    %
    function self = newton_class( varargin )
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
    function setup( self, f, Df )
      self.fun  = f;
      self.Dfun = Df;
    end
    %
    % solve the problem f(x) = 0 starting from point x0
    %
    function x = solve( self, f, Df, x0 )
      self.setup(f,Df);
      self.flag = false;
      x = x0;    % initial approximation
      for iter=1:self.max_iter
        fx  = feval( self.fun,  x ); % compute f(x);
        Dfx = feval( self.Dfun, x ); % compute f`(x);
        if Dfx == 0 % avoid division by 0
          if self.verbose == 'iter'
            fprintf( 'find Df(x) == 0\n');
          end
          break;
        end
        % perform Newton step
        x = x - fx/Dfx;
        % print iteration
        if self.verbose == 'iter'
          fprintf( 'iter=%2d x=%g f(x)=%g\n', iter, x, fx );
        end
        % check convergence
        if abs(fx) <= self.tol
          if self.verbose == 'iter'
            fprintf('convergence reached at iter %d\n',iter);
          end
          self.flag = true;
          break;
        end
      end
      self.iter   = iter;
      self.x_last = x;
    end
  end
end
