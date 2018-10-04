clear all 
clc

%declare basic variables 
gamma = 1.4;
p4 = 2.0; 
p1 = 1.0; 
rho4 = 2.0; 
rho1 = 1.0; 
imax = 41; 
xmin = 0; 
xmax = 2.; 
dx = (xmax-xmin)/(imax-1); 
x0 = (imax-1)/2;

%set up x vector 

x = 0:dx:2;

%set up state vector

u_state = zeros(3,imax+1);

%set up initial condition

for i = 1:imax+1
  u(i) = 0;
  if i == x0
    rho(1:i) = rho4;
    p(1:i) = p4;
  else
    rho(i) = rho1;
    p(i) = p1;
  end
end





u_state(1,:) = rho;
u_state(2,:) = rho.*u; 
u_state(3,:) = (p/(gamma-1)); %e




%initialize matrices

lamda_plus = zeros(3,3);
lamda_minus = zeros(3,3);
ca = zeros(3,3);
s = zeros(3,3);
ca_inverse = zeros(3,3);
s_inverse = zeros(3,3);

%loop properties 

t = 0; 
tmax = 2; 

dt = 0.02;

u_update = zeros(3,imax);

fplus = zeros(3,imax); 
fminus = zeros(3,imax);




count = 0;
while count < 18
  
%% DO FOR I + 1/2  
 
  
  for i = 2:1:imax
    ubar = 0.5*(u_state(:,i) + u_state(:,i+1));
    averageP = (ubar(3)-(0.5.*ubar(1).*(ubar(2)/ubar(1)).^2)).*(gamma-1);
    a = sqrt((averageP.*gamma)./(ubar(1)));   
    
  
    %FORM LAMDA MATRIX 
    
    l1 = ubar(2)/(ubar(1)); 
    l2 = ubar(2)/ubar(1) + a;
    l3 = ubar(2)/ubar(1) -a;
    
    l1p = 0.5*(l1+abs(l1));
    l2p = 0.5*(l2+abs(l2));
    l3p = 0.5*(l3+abs(l3));
    
    l1m = 0.5*(l1-abs(l1));
    l2m = 0.5*(l2-abs(l2));
    l3m = 0.5*(l3-abs(l3));
    
    lamda_plus(1,1) = l1p;
    lamda_plus(2,2) = l2p;
    lamda_plus(3,3) = l3p;
    
    lamda_minus(1,1) = l1m;
    lamda_minus(2,2) = l2m;
    lamda_minus(3,3) = l3m;
    
    %%FORM Ca, S and CaInv, S inv 
    
    %FOR S
    alpha = 0.5*((ubar(2)/(ubar(1))).^2);
    beta = gamma-1;
    
    s(1,1) = 1;
    s(2,1) = -ubar(2)/ubar(1)/ubar(1);
    s(2,2) = 1/ubar(1);
    s(3,1) = alpha*beta;
    s(3,2) = -(ubar(2)/ubar(1))*beta;
    s(3,3) = beta;
    
    %FOR S INVERSE 
    
     
    s_inverse(1,1) = 1;
    s_inverse(2,1) = -ubar(2)/ubar(1);
    s_inverse(2,2) = ubar(1);
    s_inverse(3,1) = alpha;
    s_inverse(3,2) = (ubar(2)/ubar(1))*ubar(1);
    s_inverse(3,3) = 1/beta;
    
    
    %FOR Ca
    
    ca(1,1) = 1; 
    ca(1,3) = -1./(a.^2);
    ca(2,2) = ubar(1)*a;
    ca(2,3) = 1; 
    ca(3,2) = -ubar(1)*a;
    ca(3,3) = 1;
    
    %FOR Ca inv 
    
    ca_inverse(1,1) = 1; 
    ca_inverse(1,2) = 1./(2*a.^2);
    ca_inverse(1,3) = 1./(2*a.^2);
    ca_inverse(2,2) = 1./(2*ubar(1)*a);
    ca_inverse(2,3) = -(1./(2*ubar(1)*a));
    ca_inverse(3,2) = 0.5;
    ca_inverse(3,3) = 0.5;
    
    aplus_p = s_inverse*ca_inverse*lamda_plus*ca*s;
    
    
    aminus_p = s_inverse*ca_inverse*lamda_minus*ca*s;
    
    fplus(:,i) = aplus_p*u_state(:,i)+aminus_p*u_state(:,i+1);
  end
  count = count + 1;
end
  
  
%     %clear variables
    
    
      
    
  

% %% FOR I MINUS HALF 
%     
%     
%     
%     ubarm = 0.5*(u_state(:,i) + u_state(:,i-1));
%     averageP = (ubarm(3)-(0.5.*ubarm(1).*(ubarm(2)/ubarm(1)).^2)).*(gamma-1);
%     a = sqrt((averageP.*gamma)./(ubarm(1)));   
%     
%     
%     
%     
%     
%   
%     %FORM LAMDA MATRIX 
%     
%     l1m = ubarm(2)/(ubarm(1)); 
%     l2m = ubarm(2)/ubarm(1) + a;
%     l3m = ubarm(2)/ubarm(1) -a;
%     
%     l1pm = 0.5*(l1m+abs(l1m));
%     l2pm = 0.5*(l2m+abs(l2m));
%     l3pm = 0.5*(l3m+abs(l3m));
%     
%     l1mm = 0.5*(l1m-abs(l1m));
%     l2mm = 0.5*(l2m-abs(l2m));
%     l3mm = 0.5*(l3m-abs(l3m));
%     
%     lamda_plusm(1,1) = l1pm;
%     lamda_plusm(2,2) = l2pm;
%     lamda_plusm(3,3) = l3pm;
%     
%     lamda_minusm(1,1) = l1mm;
%     lamda_minusm(2,2) = l2mm;
%     lamda_minusm(3,3) = l3mm;
%     
%     %%FORM Ca, S and CaInv, S inv 
%     
%     %FOR S
%     alpham = 0.5*((ubarm(2)/(ubarm(1))).^2);
%     beta = gamma-1;
%     
%     sm(1,1) = 1;
%     sm(2,1) = -ubarm(2)/ubarm(1)/ubarm(1);
%     sm(2,2) = 1/ubarm(1);
%     sm(3,1) = alpham*beta;
%     sm(3,2) = -(ubarm(2)/ubarm(1))*beta;
%     sm(3,3) = beta;
%     
%     %FOR S INVERSE 
%     
%      
%     s_inversem(1,1) = 1;
%     s_inversem(2,1) = -ubarm(2)/ubarm(1);
%     s_inversem(2,2) = ubarm(1);
%     s_inversem(3,1) = alpham;
%     s_inversem(3,2) = (ubarm(2)/ubarm(1))*ubarm(1);
%     s_inversem(3,3) = 1/beta;
%     
%     
%     %FOR Ca
%     
%     cam(1,1) = 1; 
%     cam(1,3) = -1./(a.^2);
%     cam(2,2) = ubarm(1)*a;
%     cam(2,3) = 1; 
%     cam(3,2) = -ubarm(1)*a;
%     cam(3,3) = 1;
%     
%     %FOR Ca inv 
%     
%     ca_inversem(1,1) = 1; 
%     ca_inversem(1,2) = 1./(2*a.^2);
%     ca_inversem(1,3) = 1./(2*a.^2);
%     ca_inversem(2,2) = 1./(2*ubarm(1)*a);
%     ca_inversem(2,3) = -(1./(2*ubarm(1)*a));
%     ca_inversem(3,2) = 0.5;
%     ca_inversem(3,3) = 0.5;
%     
%     aplus_m = s_inversem*ca_inversem*lamda_plusm*cam*sm;
%     aminus_m = s_inversem*ca_inversem*lamda_minusm*cam*sm;
%     
%     fminus(:,i) = aplus_m*u_state(:,i-1)+aminus_m*u_state(:,i);
% 
%     
%     
%     
%     u_update(:,i) = u_state(:,i)-(dt/dx)*(fplus(:,i)-fminus(:,i));
%     u_state(:,imax+1) = u_state(:,imax); 
% 
%     
%     t = t+ dt;
%   end
%     u_update(:,imax+1) = u_update(:,imax);
%     u_state(:,i) = u_update(:,i);
%     count = count + 1;
% end
  

% plot(x,u_update(1,1:41));
  
  
  
  








