clear all
clc

%BASIC PROPERTIES
gamma = 1.4;
xmin = 0.0;
xmax = 2.0;
imax = 41;
dx = (xmax-xmin)/(imax-1);
x0 = 1.0; %initial location diaphragmn
p1 = 1.0;
p4 = 2.0;
rho1 = 1.0;
rho4 = 4.0;
lamda1 = 0;
lamda2 = 0;
lamda3 = 0;


%VECTORS
x = 0:dx:2;
x = x.';
rho = zeros(imax,1);
p = zeros(imax,1);
u = zeros(imax,1); 
ubar = zeros(imax,1);
usolution = zeros(3,1);



%MATRICES

lamda_matrix = zeros(3,3);
s = zeros(3,3);
c_a = zeros(3,3);
x_matrix = zeros(3,3);
lamda_matrix_plus = zeros(3,3);
lamda_matrix_minus = zeros(3,3);
abar_plus = zeros(3,3);
abar_minus = zeros(3,3);
%INITIAL CONDITIONS
locateX0 = find(x==1.0);

p(1:locateX0) = p4;
p(locateX0:imax) = p1;

rho(1:locateX0) = rho4;
rho(locateX0:imax) = rho1;

%LOOP PROPERTIES
t = 0; %start time
tmax = 1; %end time

for t = 1 
  
  %calculate a
  for i = 1: imax-1
  a = sqrt(gamma.*p(i)./rho(i));
  ubar(i) = 0.5.*(u(i)+ u(i+1));
  
  lamda1 = ubar(i);
  lamda2 = ubar(i)+a;
  lamda3 = ubar(i)-a;
  
  lamda1plus = 0.5*(lamda1+abs(lamda1));
  lamda2plus = 0.5*(lamda2+abs(lamda2));
  lamda3plus = 0.5*(lamda3+abs(lamda3));
  
  
  lamda1minus = 0.5*(lamda1-abs(lamda1));
  lamda2minus = 0.5*(lamda2-abs(lamda2));
  lamda3minus = 0.5*(lamda3-abs(lamda3));
  
  %FILL LAMDA MATRICES
  lamda_matrix_plus(1,1) = lamda1plus;
  lamda_matrix_plus(2,2) = lamda2plus;
  lamda_matrix_plus(3,3) = lamda3plus;
  
  
  lamda_matrix_minus(1,1) = lamda1minus;
  lamda_matrix_minus(2,2) = lamda2minus;
  lamda_matrix_minus(3,3) = lamda3minus;
  
  
  
  
  
  alpha = (ubar(i).^2./2.);
  beta = gamma-1;
  
  %FILL S 
  s(1,1) = 1;
  s(2,1) = -u(i)./rho(i);
  s(2,2) = 1./rho(i);
  s(3,1) = alpha*beta;
  s(3,2) = -u(i)*beta;
  s(3,3) = beta;
  
  %FIND INVERSE S
  s_inverse = inv(s);
  
  %FILL CA
  ca(1,1) = 1;
  ca(1,3) = -1./(a.^2.);
  ca(2,2) = rho(i)*a;
  ca(2,3) = 1;
  ca(3,2) = -rho(i)*a;
  ca(3,3) = 1;
  
  %FIND INVERSE CA
  ca_inverse = inv(ca);
  
  %CALCULATE X 
  x_matrix = ca.*s;
  x_matrix_inverse = inv(x_matrix);
  
  
  %CALCULATE A BAR  
  
  
  abar_minus = x_matrix*lamda_matrix_minus*x_matrix_inverse;
  abar_plus = x_matrix*lamda_matrix_plus*x_matrix_inverse;
  
  
  for i = 1:imax-1
    fplus(i) = 
  
  end
  
  
  
  
  
 
  
  
  
  
  
  
  
  
  end
  
  
  
  
  

end








