function [f1, f2, f3, interval1, interval2, interval3, degree] = generate_vars()
  f1 = @(x) sin(x);
  f2 = @(x) cos(x);
  f3 = @(x) atan(x);

  interval1 = [-pi/2, pi/2];
  interval2 = [-pi/2, pi/2];
  interval3 = [-1, 1];

  degree = 5;
