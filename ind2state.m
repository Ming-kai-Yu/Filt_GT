function x = ind2state(index, base)
  % inverse conversion of state2ind
  x = zeros(4,1);
  num = index-1;
  x(1)= floor(num/base^3);
  num = num-x(1)*base^3;
  x(2) = floor(num/base^2);
  num = num-x(2)*base^2;
  x(3) = floor(num/base);
  x(4) = num-x(3)*base;
end