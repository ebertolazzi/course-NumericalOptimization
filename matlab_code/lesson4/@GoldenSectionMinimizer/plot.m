function plot( self, y0, dy )
  y = y0;
  hold on;
  for k=1:size(self.history,2)
    a = self.history(1,k);
    b = self.history(2,k);
    plot( [a,b], [y,y], '-b', 'LineWidth', 2 );
    plot( [a,a], [y-dy/3,y+dy/3], '-b', 'LineWidth', 2 );
    plot( [b,b], [y-dy/3,y+dy/3], '-b', 'LineWidth', 2 );
    y = y+dy;
  end
end
