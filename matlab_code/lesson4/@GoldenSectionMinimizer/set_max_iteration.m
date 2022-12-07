%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%
% Change the default maximum number of iteration
%
function set_max_iteration( self, max_iter )
  if length(max_iter) > 1 || ~isinteger(max_iter)
    error('GoldenSectionMinimizer:set_max_iteration, expected a scalar  integer\n');
  end
  if max_iter < 0 || max_iter > 10000
    error('GoldenSectionMinimizer:set_max_iteration, bad number of iterator %d\n',max_iter);
  end
  self.max_iter = max_iter;
end
