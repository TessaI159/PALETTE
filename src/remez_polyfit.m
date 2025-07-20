function [coeffs, err] = remez_polyfit(f, interval, degree)
  n = degree;
  xi = cheby_nodes(n + 2, interval);
  k = 0:(n+1);
  alternation = (-1).^k;
  pows = 1:degree;
  A = zeros(n+2, n+2);
  for i = 1:(n+2)
    A(i, 1) = 1;
    A(i, 2:n+1) = xi(i).^pows;
    A(i, end) = alternation(i);
  endfor
  A;
  r = f(xi(:));
  sol = A \ r;
endfunction

function [nodes] = cheby_nodes(n, interval)
  k = 0:(n-1);
  nodes_init = cos(((2*k + 1)*pi)/(2*(n)));
  a = interval(1);
  b = interval(2);
  nodes = 0.5*(a + b) + 0.5*(b - a)*nodes_init;
endfunction

function [dx] = prime(f, x, tdelta=1e-8)
  h = 0.1;
  dx = (f(x+h)- f(x))/h;
  do
    h /= 10;
    new_dx = (f(x + h) - f(x))/h;
    delta = abs(new_dx - dx);
    dx = new_dx;
  until (delta < tdelta)
endfunction

function [independent, dependent] = find_lgst_abs_extrema(f, gridsize, interval)
  grid = linspace(interval(1), interval(2), gridsize);
  old_dx = prime(f, grid(1));
  for i = 2:gridsize
    new_dx = prime(f, grid(i));
    if (new_dx * old_dx <= 0)
      
    endif
  endfor
endfunction

function [x, zero] = find_zero(f, initial, num_iter)
  if (prime(f, initial) == 0)
    x = initial;
    zero = f(x);
    return;
  endif
  i = 1;
  do
    x = x - f(x)/prime(f, x);
  until (f(x) == 0 | i == num_iter);
endfunction
