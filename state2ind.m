function index = state2ind(x, base)
   % n - total copy number of all species
   % base = n+1
   % state (0,0,0,0) maps to index 1, (0,0,0,1) to 2,
   % (0, 0, 0, n) to base
   % (0, 0, 1, 0) to base + 1, etc.
   index = x(1)*base^3 + x(2)*base^2 + x(3)*base + x(4) + 1;
end