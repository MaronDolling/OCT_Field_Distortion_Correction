function d = point_to_line_distance(pt, l1, l2)
% Calculates the distance of a point pt to a line crossing points l1 and l2
  a = l1 - l2;
  b = pt - l2;
  c = cross(a, b, 2);
  d = sqrt(sum(c.^2, 2)) ./ sqrt(sum(a.^2, 2));
end
