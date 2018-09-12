classdef MinimizationConjugateGradient < MinimizationND
  %
  % Description
  % -----------
  % Minimization of a nonlinear multiuvariate functions using nonlinbear conjugate gradient mnethods.
  % The algorithm is describbed in the references.
  %
  % References
  % ----------
  %
  % @article{2006,
  %     author   = {Hager, William W. and Zhang, Hongchao},
  %     journal  = {Pacific journal of Optimization},
  %     keywords = {bibliot-11-08-07},
  %     number   = {1},
  %     pages    = {35--58},
  %     title    = {A survey of nonlinear conjugate gradient methods},
  %     volume   = {2},
  %     year     = {2006}
  % }
  % @article{2006,
  %     author  = {Hager, William W. and Zhang, Hongchao},
  %     title   = {Algorithm 851: CG_DESCENT, a Conjugate Gradient Method with Guaranteed Descent},
  %     journal = {ACM Trans. Math. Softw.},
  %     volume  = {32},
  %     number  = {1},
  %     year    = {2006},
  %     pages   = {113--137},
  %     doi     = {10.1145/1132973.1132979},
  % }
  %
  % Author: Enrico Bertolazzi
  %

  properties (SetAccess = private, Hidden = true)
    method          % selected nonliner CG method
    method_names    % list of available methods
    angle_too_small %
  end

  methods (Hidden = true)
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = HS( ~, g1, g0, d0, ~ )
      % M.R.Hestenes and E.L.Stiefel
      % Methods of conjugate gradients for solving linear systems,
      % J. Research Nat. Bur. Standards, 49 (1952), pp. 409–436
      y    = g1 - g0;
      beta = dot(g1,y)/dot(d0,y);
      d    = -g1 + beta * d0;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = HSplus( ~, g1, g0, d0, ~ )
      % M. J. D. Powell
      % Nonconvex minimization calculations and the conjugate gradient method,
      % Numerical Analysis (Dundee, 1983), Lecture Notes in Mathematics, Vol. 1066,
      % Springer-Verlag, Berlin, 1984, pp. 122–141.
      y    = g1 - g0;
      beta = max(0,dot(g1,y)/dot(d0,y));
      d    = -g1 + beta * d0;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = FR( ~, g1, g0, d0, ~ )
      % R. Fletcher and C. Reeves
      % Function minimization by conjugate gradients,
      % Comput. J., 7 (1964), pp. 149–154.
      beta = dot(g1,g1)/dot(g0,g0);
      d    = -g1 + beta * d0;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = PRP( ~, g1, g0, d0, ~ )
      % E. Polak and G. Ribiere
      % Note sur la convergence de directions conjugees,
      % Rev. Francaise Informat Recherche Opertionelle, 3e Annee 16 (1969), pp. 35–43.
      %
      % B. T. Polyak,
      % The conjugate gradient method in extreme problems,
      % USSR Comp. Math. Math. Phys., 9 (1969), pp. 94–112.
      y    = g1 - g0;
      beta = dot(g1,y)/dot(g0,g0);
      d    = -g1 + beta * d0;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = PRPplus( ~, g1, g0, d0, ~ )
      % J. C. Gilbert and J. Nocedal
      % Global convergence properties of conjugate gradient methods for optimization,
      % SIAM J. Optim., 2 (1992), pp. 21–42.
      y    = g1 - g0;
      beta = max(0,dot(g1,y)/dot(g0,g0));
      d    = -g1 + beta * d0;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = aPRP( ~, g1, g0, ~, s )
      % N. Andrei
      % Another nonlinear conjugate gradient algorithm with sufficient descent conditions for unconstrained optimization
      % ICI Technical Reports, November 22, 2006
      y    = g1-g0;
      beta = ( dot(y,g1) - dot(y,y)*dot(s,g1)/dot(g0,g0) )/dot(y,s); % ORIGINAL
      d    = -g1 + beta * s;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = CCOMB( ~, g1, g0, ~, s )
      % N. Andrei
      % New hybrid conjugate gradient algorithms as a convex combination of PRP and DY for unconstrained optimization
      % ICI Technical Report, October 1, 2007.
      % CHANGED s <- d0
      % check restart criteria of Powell
      y     = g1 - g0;
      ys    = dot(y,s);
      yg1   = dot(y,g1);
      g00   = dot(g0,g0);
      g11   = dot(g1,g1);
      prp   = yg1/g00;
      dy    = g11/ys;
      theta = max( 0, min( 1, yg1*(ys-g00)/(yg1*ys-g00*g11) ) );
      beta  = (1-theta)*prp + theta * dy;
      d     = -g1 + beta * s;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = NDOMB( ~, g1, g0, ~, s )
      % N. Andrei
      % New hybrid conjugate gradient algorithms for unconstrained optimization
      % Encyclopedia of Optimization, 2nd Edition, C.A. Floudas, and P. Pardalos (Eds.), August 2008, Entry 761.
      y     = g1 - g0;
      ys    = dot(y,s);
      yg1   = dot(y,g1);
      g00   = dot(g0,g0);
      g11 = dot(g1,g1);
      prp   = yg1/g00;
      dy    = g11/ys;
      theta = max( 0, min( 1, (yg1*ys - g00*(yg1-dot(s,g1)))/(yg1*ys-g00*g11) ) );
      beta  = (1-theta)*prp + theta * dy;
      d     = -g1 + beta * s;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = CD( ~, g1, g0, d0, ~ )
      % R. Fletcher
      % Practical Methods of Optimization vol. 1: Unconstrained Optimization,
      % John Wiley & Sons, New York, 1987.
      beta = -dot(g1,g1)/dot(g0,d0);
      d    = -g1 + beta * d0;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = LS( ~, g1, g0, d0, ~ )
      % Y. Liu and C. Storey
      % Efficient generalized conjugate gradient algorithms, Part 1: Theory,
      % J. Optim. Theory Appl., 69 (1991), pp. 129–137.
      y    = g1 - g0;
      beta = -dot(g1,y)/dot(g0,d0);
      d    = -g1 + beta * d0;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = DY( ~, g1, g0, d0, ~ )
      % Y. H. Dai and Y. Yuan
      % A nonlinear conjugate gradient method with a strong global convergence property
      % SIAM J. Optim., 10 (1999), pp. 177–182.
      y    = g1 - g0;
      beta = dot(g1,g1)/dot(y,d0);
      d    = -g1 + beta * d0;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = DL( ~, g1, g0, ~, s )
      % Y.H. Dai and L.Z. Liao
      % New conjugacy conditions and related nonlinear conjugate gradient methods
      % Appl. Math. Optim., 43 (2001), pp. 87-101
      y    = g1 - g0;
      t    = 0.1;
      beta = dot(g1,y-t*s)/dot(y,s);
      d    = -g1 + beta * s;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = DLplus( ~, g1, g0, ~, s )
      % H. Yabe and M. Takano
      % Global convergence properties of nonlinear conjugate gradient methods with modified secant conditions
      % Computational Optimization and Applications, 28 (2004), pp.203-225.
      y    = g1 - g0;
      t    = 0.1;
      ys   = dot(y,s);
      beta = max(0,dot(g1,y)/ys) - t*dot(g1,s)/ys;
      d    = -g1 + beta * s;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = ACGA( ~, g1, g0, ~, s )
      % N. Andrei
      % Another nonlinear conjugate gradient algorithm for unconstrained optimization
      % ICI Technical Report, May 16, 2007.
      y    = g1 - g0;
      ys   = dot(y,s);
      beta = (dot(y,g1)/ys) * ( 1 - dot(s,g1)/ys ); % changed s <= d0
      d    = -g1 + beta * s;
      % check if necessary to reset direction
      if dot(d,g1) > -1e-3 * norm(d) * norm(g1)
        d = -g1;
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = ACGAplus( ~, g1, g0, ~, s )
      % N. Andrei
      % Another nonlinear conjugate gradient algorithm for unconstrained optimization
      % ICI Technical Report, May 16, 2007.
      y    = g1 - g0;
      ys   = dot(y,s);
      beta = max(0,dot(y,g1)/ys) * ( 1 - dot(s,g1)/ys );
      d    = -g1 + beta * s;
      % check if necessary to reset direction
      if dot(d,g1) > -1e-3 * norm(d) * norm(g1)
        d = -g1;
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = CGSD( ~, g1, g0, ~, s )
      % N. Andrei, A Dai-Yuan
      % conjugate gradient algorithm with sufficient descent and conjugacy conditions for unconstrained optimization.
      % Applied Mathematics Letters, vol.21, 2008, pp.165-171.
      y     = g1-g0;
      ys    = dot(y,s);
      yg1   = dot(y,g1);
      g11   = dot(g1,g1);
      theta = g11/yg1;
      beta  = (g11-yg1*dot(s,g1)/ys)/ys;
      d     = -theta * g1 + beta * s;
      % check if necessary to reset direction
      if dot(d,g1) > -1e-3 * norm(d) * sqrt(g11)
        d = -g1;
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = BM( ~, g1, g0, ~, s )
      % E. Birgin and J.M. Martínez
      % A spectral conjugate gradient method for unconstrained optimization
      % Applied Mathematics and Optimization, 43, pp.117-128, 2001.
      y     = g1-g0;
      ys    = dot(y,s);
      theta = dot(s,s)/ys;
      beta  = dot(theta*y-s,g1)/ys;
      d     = -theta * g1 + beta * s;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = BMplus( ~, g1, g0, ~, s )
      % E. Birgin and J.M. Martínez
      % A spectral conjugate gradient method for unconstrained optimization
      % Applied Mathematics and Optimization, 43, pp.117-128, 2001.
      y     = g1-g0;
      ys    = dot(y,s);
      theta = dot(s,s)/ys;
      beta  = (max(0,dot(theta*y,g1))-dot(s,g1))/ys;
      d     = -theta * g1 + beta * s;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = N( ~, g1, g0, d0, s )
      % W. W. Hager and H. Zhang
      % A new conjugate gradient method with guaranteed descent and an efficient line search
      % SIAM J. Optim., 16(1), 170–192. 2003
      y    = g1 - g0;
      z    = y - (2*dot(y,y)/dot(d0,y))*d0;
      eta  = -1/(norm(s)*min(0.01,norm(g0)));
      beta = max( eta, dot(z,g1)/dot(s,y) );
      d    = -g1 + beta * s;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = HyTAS( ~, g1, g0, d0, ~ )
      % D. Touati-Ahmed and C. Storey
      % Efficient hybrid conjugate gradient techniques,
      % J. Optim. Theory Appl., 64 (1990), pp. 379–397.
      g0g0    = dot(g0,g0);
      betaPRP = dot(g1,g1-g0)/g0g0;
      betaFR  = dot(g1,g1)/g0g0;
      if 0 <= betaPRP && betaPRP <= betaFR
        beta = betaPRP;
      else
        beta = betaFR;
      end
      d = -g1 + beta * d0;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = HyHS( ~, g1, g0, d0, ~ )
      % Y. F. Hu and C. Storey
      % Global convergence result for conjugate gradient methods,
      % J. Optim. Theory Appl., 71 (1991), pp. 399–405.
      g0g0    = dot(g0,g0);
      betaPRP = dot(g1,g1-g0)/g0g0;
      betaFR  = dot(g1,g1)/g0g0;
      beta    = max(0,min(betaPRP,betaFR));
      d       = -g1 + beta * d0;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = HyNG( ~, g1, g0, d0, ~ )
      % J. C. Gilbert and J. Nocedal
      % Global convergence properties of conjugate gradient methods for optimization,
      % SIAM J. Optim., 2 (1992), pp. 21–42.
      g0g0    = dot(g0,g0);
      betaPRP = dot(g1,g1-g0)/g0g0;
      betaFR  = dot(g1,g1)/g0g0;
      beta    = max(-betaFR,min(betaPRP,betaFR));
      d       = -g1 + beta * d0;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = HyDY( ~, g1, g0, d0, ~ )
      % Y. H. Dai and Y. Yuan
      % An efficient hybrid conjugate gradient method for unconstrained optimization,
      % Ann. Oper. Res., 103 (2001), pp. 33–47.
      beta = dot(g1,g1)/max( dot(d0,g1-g0), -dot(g0,d0));
      d    = -g1 + beta * d0;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function d = dirEval( self, g1, g0, d0, s )
      % Compute the value of `beta` for the selected nonlinear CG method
      switch ( self.method )
      case 'GRAD';  d = -g1; % standard gradient method
      case 'ACGA';  d = self.ACGA( g1, g0, d0, s );
      case 'ACGA+'; d = self.ACGAplus( g1, g0, d0, s );
      case 'BM';    d = self.BM( g1, g0, d0, s );
      case 'BM+';   d = self.BMplus( g1, g0, d0, s );
      case 'CCOMB'; d = self.CCOMB( g1, g0, d0, s );
      case 'CD';    d = self.CD( g1, g0, d0, s );
      case 'CGSD';  d = self.CGSD( g1, g0, d0, s );
      case 'DL';    d = self.DL( g1, g0, d0, s );
      case 'DL+';   d = self.DLplus( g1, g0, d0, s );
      case 'DY';    d = self.DY( g1, g0, d0, s );
      case 'FR';    d = self.FR( g1, g0, d0, s );
      case 'HS';    d = self.HS( g1, g0, d0, s );
      case 'HS+';   d = self.HSplus( g1, g0, d0, s );
      case 'HyDY';  d = self.HyDY( g1, g0, d0, s );
      case 'HyHS';  d = self.HyHS( g1, g0, d0, s );
      case 'HyNG';  d = self.HyNG( g1, g0, d0, s );
      case 'HyTAS'; d = self.HyTAS( g1, g0, d0, s );
      case 'LS';    d = self.LS( g1, g0, d0, s );
      case 'N';     d = self.N( g1, g0, d0, s );
      case 'NDOMB'; d = self.NDOMB( g1, g0, d0, s );
      case 'aPRP';  d = self.aPRP( g1, g0, d0, s );
      case 'PRP';   d = self.PRP( g1, g0, d0, s );
      case 'PRP+';  d = self.PRPplus( g1, g0, d0, s );
      otherwise
        error('MinimizationCG, dirEval no method selected (this should not be happen)');
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end

  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = MinimizationConjugateGradient( fun, ls )
      %
      % fun = function to be minimized
      % ls  = linesearch used, can be LinesearchArmijo, LinesearchWolfe, LinesearchGoldenSection
      self@MinimizationND( fun, ls );
      self.method_names = { ...
        'GRAD', ...
        'ACGA','ACGA+', ...
        'BM', 'BM+', ...
        'CCOMB', 'CD', 'CGSD', ...
        'DL', 'DL+', 'DY', ...
        'FR', ...
        'HS','HS+', ...
        'HyDY', 'HyHS', 'HyNG', 'HyTAS', ...
        'LS', ...
        'N', 'NDOMB', ...
        'aPRP', 'PRP','PRP+' ...
      };
      self.method          = self.method_names{1};
      self.angle_too_small = cos(pi/2+pi/180); % 90-1 degree
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function n = numOfMethods( self )
      % return the number of methods available
      n = length(self.method_names);
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function name = activeMethod( self )
      % return active CG method used in minimization
      name = self.method;
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function selectByNumber( self, k )
      % select the CG method by number
      if k < 1 || k > length(self.method_names)
        error('MinimizationCG, selectByNumber, k=%d out of range',k);
      end
      self.method = self.method_names{k};
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function selectByName( self, name )
      % select the CG method by its name
      if ischar(name)
        for k=1:length(self.method_names)
          if strcmp(name,self.method_names{k})
            self.method = name;
            return;
          end
        end
        error('MinimizationCG, selectByName, name=%s not found',name);
      else
        error('MinimizationCG, selectByName, expected string as arument, found %s',class(name));
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function [x,converged] = minimize( self, x0 )
      % generic onjugate gradient minimization algorithm
      % x0 = initial point
      %
      self.n_fail = 0;
      x           = x0;
      alpha       = 1;
      converged   = false;
      if self.save_iterate
        self.x_history = x(:);
      end
      %
      for iter=1:self.max_iter
        self.iter = iter;
        %
        % gradient of the function
        g1 = self.funND.grad( x ).';
        %
        % check if converged
        norm_inf_g1 = norm(g1,inf);
        converged = norm_inf_g1 < self.tol;
        if converged
          if self.verbose
            fprintf('solution found, ||grad f||_inf = %g < %g (tolerance)\n', ...
                     norm_inf_g1, self.tol );
          end
          break;
        end
        %
        % only for debug
        if self.verbose
          fprintf('[%s] iter = %5d ||grad f||_inf = %12.6g', ...
                  self.method, self.iter, norm_inf_g1 );
        end
        %
        % build search direction
        %
        if self.iter == 1
          d = -g1; % first iteration, search direction is -gradient
        elseif abs(dot(g0,g1)) >= 0.2 * dot(g1,g1) % check restart criteria of Powell
          if self.verbose
            fprintf(2,' [Powell restart] ' );
          end
          d = -g1; % reset direction
        else
          d = self.dirEval( g1, g0, d, s ); % s = x1 - x0 last step!
          % check if the direction is descending,
          % if direction is >= 89 degree from gradient reset to gradient directions
          % use >= to catch d == 0
          norm_2_d  = norm(d);
          norm_2_g1 = norm(g1);
          if dot( d, -g1 ) <= norm_2_d * norm_2_g1 * self.angle_too_small
            if self.verbose
              fprintf(2,' [angle ~90, reset] ' );
            end
            d = -g1; % reset direction
          end
        end
        % last check of direction search
        if norm(d,inf) == 0
          error('MinimizationCG, bad direction d == 0\n');
        end
        %
        % minimize along search direction
        [ x, alpha, ok ] = self.step1D( x, d, 10*alpha );
        %if ~ok
        %  % step failed try to use gradient direction
        %  d = -g1;
        %  [x,alpha,ok] = self.step1D( x, d, alpha );
          if ~ok
            % cannot advance see if accept a low precision solution
            fprintf(2,'\nMinimizationCG, step1D failed\n');
            return;
          end
        %end
        %
        if self.verbose
          fprintf(' alpha = %8.4g\n', alpha);
        end
        %
        % save old gradient and old step
        g0 = g1;
        s  = alpha*d;
      end
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
