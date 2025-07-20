function [x, zero] = find_zero(f, initial, num_iter)
  x = initial;
  if (prime(f, initial) == 0)
    zero = f(x);
    return;
  endif
  i = 1;
  do
    printf("prime(f, %.8f): %.8f\n", x, prime(f, x));
    x = x - f(x)/prime(f, x);
    i++;
    printf("x: %.5f, f(x): %.5f\n", x, f(x));
  until (f(x) == 0 | i == num_iter);
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
