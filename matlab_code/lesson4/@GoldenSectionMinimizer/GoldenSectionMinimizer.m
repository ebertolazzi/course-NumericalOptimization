classdef GoldenSectionMinimizer < handle

  properties (SetAccess = private, Hidden = true)
    tau      % golden section ratio (sqrt(5)-1)/2
    tol      % tolerance for |b-a|
    max_iter % maximum number of admitted iterate for Golde search
    history  % 2xn matrix with the hystory of the computed intervals
  end

  methods
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    % Construcot method called when object is instantiated
    %
    % GS = GoldenSectionMinimizer();
    %
    function self = GoldenSectionMinimizer()
      self.tau      = (sqrt(5)-1)/2;
      self.tol      = 1e-4;
      self.max_iter = 8;
      self.history  = [];
    end
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    plot( self, y0, dy );
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    set_max_iteration( self, max_iter );
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    set_tolerance( self, tol )
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    [a,b] = minimize( self, fun, a_in, b_in );
    %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end
end
