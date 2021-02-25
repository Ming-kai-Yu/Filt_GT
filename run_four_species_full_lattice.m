%% building the CME matrix

tic;

x0 = [0; 0; 0; 20];
c = [1; 1.5; 1.2; 1.5];

nu1 = [-1; 1; 0; 0];
nu2 = [0; -1; 1; 0];
nu3 = [0; 0; -1; 1];
nu4 = [1; 0; 0; -1];

nu = [nu1, nu2, nu3, nu4];


base = sum(x0)+1;
num_node = base^4;

%global A_global
A = sparse(num_node, num_node);
ind_state = zeros(num_node, 4);

for i=1:num_node
    x = ind2state(i,base);
    A(i,i) = -sum(prop(x, c));
    ind_state(i,:)=x';
    for reac=1:4
       x_in = x - nu(:,reac);
         if prod (x_in>=0 & x_in<= sum(x0))
             j = state2ind(x_in, base);
             props = prop(x_in, c);
             A(i,j)=props(reac);
         end
    end
end

%% run ODE
T = 3
tspan = [0, T];
index0 = state2ind(x0, base);
p0 = zeros(num_node, 1);
p0(index0) = 1;

[t, p] = ode23(@(t, p) four_species_cme_full_lattice(t, p, A), tspan, p0);
%%
p_final = p(end,:)';
[p1, p2, p3, p4] = joint2margnl(p_final, base);
p_margnl = [p1, p2, p3, p4]
% It is tricky to find the *stationary distribution* by solving Ax=0.
% I would expect null(A)=sum(x0)+1, and numericall it seems true.
% I would expect each vector in the null space corresponds to 
% one invariant space that x1+x2+x3+x4=k, where k=0,1,...,sum(x0).
% However, eigen vectors in null space doesn't seem easy to be aligned
% to the invariant space.
p34 = get_jointx3x4(p_final, base);
surf(p34)

%% filtering
x3 = 2; x4 = 3;
[pi1, pi2] = p1p2_given_x3x4(p_final, base, x3, x4);
pi_margnl = [pi1, pi2]

toc;
%% functions deal with conversion
function index = state2ind(x, base)
   % n - total copy number of all species
   % base = n+1
   % state (0,0,0,0) maps to index 1, (0,0,0,1) to 2,
   % (0, 0, 0, n) to base
   % (0, 0, 1, 0) to base + 1, etc.
   index = x(1)*base^3 + x(2)*base^2 + x(3)*base + x(4) + 1;
end

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

function a = prop(x,c)
  a = [c(1)*x(1); c(2)*x(2); c(3)*x(3); c(4)*x(4)];
end

function [p1, p2, p3, p4] = joint2margnl(p, base)
  p1 = zeros(base, 1);
  p2 = zeros(base, 1);
  p3 = zeros(base, 1);
  p4 = zeros(base, 1);
  states = 0:base-1;
  num_node = length(p);
  for i=1:num_node
      x = ind2state(i, base);
      p1(x(1)+1)=p1(x(1)+1) + p(i);
      p2(x(2)+1)=p2(x(2)+1) + p(i);
      p3(x(3)+1)=p3(x(3)+1) + p(i);
      p4(x(4)+1)=p4(x(4)+1) + p(i);
  end
end


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

function p34 = get_jointx3x4(p, base)
  p34 = zeros(base, base);
  %states = 0:base-1;
  num_node = length(p);
  for i=1:num_node
      x = ind2state(i, base);
      p34(x(3)+1, x(4)+1) = p34(x(3)+1, x(4)+1) + p(i);
  end
end