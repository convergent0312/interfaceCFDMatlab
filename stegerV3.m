clear all 
clc 
gamma = 1.4;
xmin = 0.0; 
xmax = 2.0; 
imax = 41; 
dx = (xmax-xmin)/(imax-1);
alpha = zeros(1,imax);
lamdaP = zeros(3,3);
lamdaM = zeros(3,3);
U_state = zeros(3,imax);
ubar = zeros(imax,1);
u = zeros(imax,1);
abarM = zeros(3,3);
abarP = zeros(3,3);
L1 = 0.0;
L2 = 0.0; 
L3 = 0.0;

rho4 = 2.0;
rho1 = 1.0;
p4 = 2.0;
p1 = 1.0;

p(1:20) = p4;
p(20:imax) = p1;

rho(1:20) = rho4;
rho(20:imax) = rho1;



% E = p/(gamma+1)+0.5*rho.*u.^2;



%STATE VECTORS COMPONENTS

U_state(1,1:20) = rho4;
U_state(1,20:imax) = rho1;

U_state(2,:) = 0;
U_state(3,1:20) = ((p4)/((gamma-1)*(rho4))).*(rho4);
U_state(3,20:imax) = ((p1)/((gamma-1)*(rho1))).*(rho1);


%LOOP PROPERTIES
t = 0; 
t_end = 1; 
dt = 0.1;




for t = 1:10
  
  
  for i = 2:imax-1
    a = sqrt(gamma.*p(i)./rho(i));
    ubar(i) = 0.5.*(u(i)+u(i+1));
    
    %define lamdas
    L1 = ubar(i);
    L2 = ubar(i)+a;
    L3 = ubar(i)-a;
    
    L1P = 0.5*(L1+abs(L1)); 
    L2P = 0.5*(L2+abs(L2));
    L3P = 0.5*(L3+abs(L3));
    
    L1N = 0.5*(L1-abs(L1));
    L2N = 0.5*(L2-abs(L2));
    L3N = 0.5*(L3-abs(L3));
    
    
    
    %FORM LAMDA MATRIX
    
    lamdaP(1,1) = L1P;
    lamdaP(2,2) = L2P;
    lamdaP(3,3) = L3P;
    
    lamdaM(1,1) = L1N;
    lamdaM(2,2) = L2N;
    lamdaM(3,3) = L3N;
    
    
    
    %HARD CODE THE X MATRIX    


    alpha = rho(i)./(a*sqrt(2));
    x_matrix(1,1) = 1;
    x_matrix(1,2) = alpha;
    x_matrix(1,3) = alpha;
    x_matrix(2,1) = ubar(i);
    x_matrix(2,2) = alpha.*(ubar(i)+a);
    x_matrix(2,3) = alpha.*(ubar(i)-a);
    x_matrix(3,1) = 0.5*ubar(i).^2;
    x_matrix(3,2) = alpha.*((0.5.*ubar(i).^2)+(ubar(i).*a)+(a.^2./(gamma-1)));
    x_matrix(3,3) = alpha.*((0.5.*ubar(i).^2)-(ubar(i).*a)+(a.^2./(gamma-1)));

   %FIND X INVERSE
   x_matrix_inverse = inv(x_matrix);

   %FIND ABAR
   abarP = x_matrix*lamdaP*x_matrix_inverse;
   abarM = x_matrix*lamdaM*x_matrix_inverse;

   %FIND FPLUS AND FMINUS
   fplus = abarP*U_state(:,i)+ abarM*U_state(:,i+1);

   fminus = abarP*U_state(:,i-1)+abarM*U_state(:,i);

     
  
  end
  
  
  

  
  
  
end

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  