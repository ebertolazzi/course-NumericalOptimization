%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Change the default tolerance
%
function set_tolerance( self, tol )
  if tol <= 0
    error('GoldenSectionMinimizer:set_tolerance, bad tolerance %g\n',tol);
  end
  self.tol = tol;
end
