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
rho4 = 2.0;
lamda1p = 0;
lamda2p = 0;
lamda3p = 0;


%VECTORS
x = 0:dx:2;
x = x.';
U_state = zeros(3,imax); 
U_state_update = zeros(3,imax);
ubarplus = zeros(imax,1);
u = zeros(imax,1);



s = zeros(3,3);
ca = zeros(3,3);


%MATRICES

lamda_matrix = zeros(3,3);
s = zeros(3,3);
ca = zeros(3,3);
ca_inverse = zeros(3,3);
abarplus = zeros(3,3);
abarminus = zeros(3,3);



%INITIAL CONDITIONS
locateX0 = find(x==1.0);

p(1:locateX0) = p4;
p(locateX0:imax) = p1;

rho(1:locateX0) = rho4;
rho(locateX0:imax) = rho1;


%FILL IN ELEMENTS FOR U_STATE (INITIAL)
% RHO, RHOU and E

U_state(1,:) = rho.';

U_state(2,:) = rho*u;

E = p/(gamma-1)+0.5*rho.*u.^2;

U_state(3,:) = rho*E;



%LOOP PROPERTIES

t = 0; %start time
tmax = 0.6846; %end time

% for t = 0:dt:18*dt
while t < tmax
  
  %calculate a
  for i = 2: imax-1
  a = sqrt(gamma.*p(i)./rho(i));
  dt = (0.9*dx)/a;
  ubarplus(i) = 0.5.*(u(i)+ u(i+1));  
  
  
  
  lamda1p = ubarplus(i);
  lamda2p = ubarplus(i)+a;
  lamda3p = ubarplus(i)-a;
  
   
  lamda1plus = 0.5*(lamda1p+abs(lamda1p));
  lamda2plus = 0.5*(lamda2p+abs(lamda2p));
  lamda3plus = 0.5*(lamda3p+abs(lamda3p));
  
  lamda1minus = 0.5*(lamda1p-abs(lamda1p));
  lamda2minus = 0.5*(lamda2p-abs(lamda2p));
  lamda3minus = 0.5*(lamda3p-abs(lamda3p));
  
  
  
  %FILL LAMDA MATRICES
  %PLUS
  lamda_plus(1,1) = lamda1plus;
  lamda_plus(2,2) = lamda2plus;
  lamda_plus(3,3) = lamda3plus;
  
  
  lamda_minus(1,1) = lamda1minus;
  lamda_minus(2,2) = lamda2minus;
  lamda_minus(3,3) = lamda3minus;
  
%% HARD CODE X MATRIX 

  %X AND X INVERSE FOLLOW BOOK
  alpha = rho(i)./(a*sqrt(2));

  x_matrix(1,1) = 1;
  x_matrix(1,2) = alpha;
  x_matrix(1,3) = alpha;
  x_matrix(2,1) = ubarplus(i);
  x_matrix(2,2) = alpha.*(ubarplus(i)+a);
  x_matrix(2,3) = alpha.*(ubarplus(i)-a);
  x_matrix(3,1) = 0.5*ubarplus(i).^2;
  x_matrix(3,2) = alpha.*((0.5.*ubarplus(i).^2)+(ubarplus(i).*a)+(a.^2/(gamma-1)));
  x_matrix(3,3) = alpha.*((0.5.*ubarplus(i).^2)-(ubarplus(i).*a)+(a.^2/(gamma-1)));
  
  

  beta = 1/(rho(i)*a*sqrt(2));
  
  x_matrix_inverse(1,1) = 1 - ((ubarplus(i).^2/2)*((gamma-1)/(a.^2)));
  x_matrix_inverse(1,2) = (gamma-1).*(ubarplus(i)/(a.^2));
  x_matrix_inverse(1,3) = -(gamma-1)*(1/a.^2);
  x_matrix_inverse(2,1) = beta.*( ((gamma-1).*(ubarplus(i).^2/2)) - (ubarplus(i).*a) );
  x_matrix_inverse(2,2) = beta.*(a- (gamma-1).*ubarplus(i));
  x_matrix_inverse(2,3) = beta.*(gamma-1);
  x_matrix_inverse(3,1) = beta.*( ((gamma-1).*(ubarplus(i).^2/2)) + (ubarplus(i).*a) );
  x_matrix_inverse(3,2) = beta.*(a+ (gamma-1).*ubarplus(i));
  x_matrix_inverse(3,3) = beta.*(gamma-1);
  
  
  






  %% CALCULATE A BAR  
  
 
  
  abarplus = x_matrix*lamda_plus*x_matrix_inverse;
  abarminus = x_matrix*lamda_minus*x_matrix_inverse;


  fplus = abarplus*U_state(:,i)+ abarminus*U_state(:,i+1);
  fminus = abarplus*U_state(:,i)+ abarminus*U_state(:,i+1);



  

  %FOLLOW BOOK, FLUX VECTORS EPLUS AND EMINUS 
  
%   eplus(1,1) = (2*gamma*u)+a-u; 
%   eplus(


  
  %% FDE
  
  for j = 2:imax-1
  
    U_state_update(1,j) = U_state(1,j) - (dt/dx).*(fplus-fminus);
    U_state(1,j) = U_state_update(1,j);
    
   
  end 
    

  t = t + dt; 
  
%   U_state(1,imax) = U_state(1,imax-1);
%   U_state(2,imax) = U_state(2,imax-1);
%   U_state(3,imax) = U_state(3,imax-1);
%   
%   U_state(1,1) = U_state(1,2);
%   U_state(2,1) = U_state(2,2);
%   U_state(3,1) = U_state(3,2);
  
  
 
  
  end   
  

end

density = U_state_update(1,:);
% plot(x,density);
% pause(0.1);






