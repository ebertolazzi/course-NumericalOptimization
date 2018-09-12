classdef MinimizationND < handle

  properties ( SetAccess = protected, Hidden = true )
    funND       % class function with expected methods
                % funND.eval(x)
                % funND.grad(x)
                % funND.hessian(x) (UNUSED)
    linesearch  % quasi minimizzation with expected method
                % linesearch.setFunction(fun)
                % [alpha,ok] = linesearch.search(alpha_guess)
    tol         % tolleranza per |grad|
    FD_D        % use finite 1D difference or gradient for directional derivative computation
    x_history   % list of iteration
    verbose     % if true do some additional printing and store iterations
    save_iterate
    max_iter    % massimo numero iterate ammesse
    iter
  end

  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = MinimizationND( fun, ls )
      self.tol          = 1e-3;
      self.max_iter     = 100;
      self.verbose      = true;
      self.save_iterate = false;
      self.FD_D         = false;
      self.funND        = fun;
      self.linesearch   = ls;
      self.iter         = 0;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setFunction( self, f )
      self.funND = f;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setLinesearch( self, ls )
      self.linesearch = ls;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setTolerance( self, tol )
      if tol <= 0
        error('MinimizationND, bad tolerance %g\n',tol);
      end
      self.tol = tol;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function setMaxIteration( self, max_iter )
      if length(max_iter) > 1 || ~isinteger(max_iter)
        error('MinimizationND, expected a scalar  integer\n');
      end
      if max_iter < 0 || max_iter > 1000000
        error('MinimizationND, bad number of iterator %d\n',max_iter);
      end
      self.max_iter = max_iter;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function iter = getIteration( self )
      iter = self.iter;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function verbose_on( self )
      self.verbose = true;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function verbose_off( self )
      self.verbose = false;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function save_iterate_on( self )
      self.save_iterate = true;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function save_iterate_off( self )
      self.save_iterate = false;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function use_FD_D( self )
      self.FD_D = true;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function no_FD_D( self )
      self.FD_D = false;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [ x1, alpha, ok ] = step1D( self, x0, d, alpha_guess )
      %
      % Perform linesearch of ND-function starting from x0 along direction d
      % grad0 is the gradient of the function at x0
      % alpha_guess is a guess for the linesearch
      %
      if norm(d,inf) == 0
        error('MinimizationND, bad direction d == 0\n');
      else
        norm_2_d = norm(d);
      end
      %
      % build the 1D function along the search direction
      fcut = Function1Dcut( self.funND, x0, d ./ norm_2_d );

      % set analitic gradient if necessary
      if self.FD_D
        fcut.use_FD_D();
      else
        fcut.no_FD_D();
      end

      %
      % do a 1D minimization
      % search an interval for minimization
      self.linesearch.setFunction( fcut );
      [ alpha, ok ] = self.linesearch.search( alpha_guess * norm_2_d );
      alpha = alpha / norm_2_d;
      %
      % check error
      if ~ok
        x1 = x0;
      else
        %
        % advance
        x1 = x0 + alpha * d;
        %
        if self.save_iterate
          self.x_history = [ self.x_history x1(:) ];
        end
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [xmin,ymin,xmax,ymax] = iterRange( self )
      if size(self.x_history,1) == 2
        xmin = min(self.x_history(1,:));
        xmax = max(self.x_history(1,:));
        ymin = min(self.x_history(2,:));
        ymax = max(self.x_history(2,:));
        dx   = xmax-xmin; dy   = ymax-ymin;
        xmin = xmin - dx; ymin = ymin - dy;
        xmax = xmax + dx; ymax = ymax + dy;
      else
        xmin = 0;
        xmax = 0;
        ymin = 0;
        ymax = 0;
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plotIter( self, last )
      if size(self.x_history,1) == 2
        hold on;
        xo = self.x_history(:,1);
        plot(xo(1),xo(2),'or');
        for k=2:size(self.x_history,2)
          xn = self.x_history(:,k);
          plot([xo(1);xn(1)],[xo(2);xn(2)],'-b');
          plot(xo(1),xo(2),'or');
          xo = xn;
        end
        xo = self.x_history(:,end);
        plot( xo(1),xo(2),'o', ...
              'MarkerSize',10,...
              'MarkerEdgeColor','b',...
              'MarkerFaceColor','black');
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function plotResidual( self, varargin )
      N = size(self.x_history,2);
      if N > 0
        r = zeros(N,1);
        for k=1:N
          r(k) = norm(self.funND.grad(self.x_history(:,k)),inf);
        end
        semilogy( 1:N, r, varargin{:} );
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
