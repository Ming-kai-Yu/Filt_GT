function [p1, p2] = p1p2_given_x3x4(p, base, x3, x4)
  p1 = zeros(base, 1);
  p2 = zeros(base, 1);
  %states = 0:base-1;
  num_node = length(p);
  for i=1:num_node
      x = ind2state(i, base);
      if (x(3)== x3 && x(4)==x4)
        p1(x(1)+1)=p1(x(1)+1) + p(i);
        p2(x(2)+1)=p2(x(2)+1) + p(i);
      end
  end
  p1 = p1/sum(p1);
  p2 = p2/sum(p2);
end