

clear all
clc

%BASIC PROPERTIES
gamma = 1.4;
xmin = 0.0;
xmax = 2.0;
imax = 41;
dx = (xmax-xmin)/(imax-1);
x0 = 1.0; %initial location diaphragmn
lamda1plus_ubarplus= 0;
lamda2plus_ubarplus= 0;

%VECTORS
x = 0:dx:2;
x = x.';
U_state = zeros(3,imax); 
U_state_update = zeros(3,imax);
u = zeros(imax,1);
ubarplus = zeros(imax,1);
ubarminus = zeros(imax,1);







%MATRICES

x_matrix_uplus = zeros(3,3);
x_matrix_inverse_uplus = zeros(3,3);

x_matrix_uminus = zeros(3,3);
x_matrix_inverse_uminus = zeros(3,3);


lamda_plus_u_plus = zeros(3,3);
lamda_minus_u_plus = zeros(3,3);
lamda_plus_u_minus = zeros(3,3);
lamda_minus_u_minus = zeros(3,3);

abarplus_uplus = zeros(3,3);
abarminus_uplus = zeros(3,3);

abarplus_uminus = zeros(3,3);
abarminus_uminus = zeros(3,3);

%INITIAL CONDITIONS
locateX0 = find(x==1.0);


%Fill in initial value for P, and RHO 
for i =  1:imax                 
    if i<=imax/2
        rho(i) = 2.0;    %density
        p(i) = 2.0;      %pressure
        u(i) = 0.0;      %velocity
    else
        rho(i) = 1.0;
        p(i) = 1.0;
        u(i) = 0.0;       
    end
end


% U is intially zero

%Calculate energy 


E = p/(gamma-1)+0.5*rho.*u.^2;



%FILL IN ELEMENTS FOR U_STATE (INITIAL)
% RHO, RHOU and RHOE

% U_state(1,:) = rho;
% U_state(2,:) = rho*u;
% U_state(3,:) = rho.*E;

u2 = rho.*u;
u3 = rho.*E;

%% START THE LOOP


%LOOP PROPERTIES
t = 0; %start time
tmax = 1; %end time
courant  = 0.9;
dt = 0.1;
tlast = 2.;

% while t< tlast

for t = 1:10
     
   
   
  
  for i = 2: imax-1
      
%% LAMDA

    a = sqrt(gamma*p(i)./rho(i));
    ubarplus(i)= (u(i) + u(i+1))/2.; %for Ubar i + 1/2
    ubarminus(i)= (u(i) - u(i-1))/2.; % for Ubar i-1/2 
    
    
    
    %SPLIT LAMDAS INTO PLUS AND MINUS
    % FOR U BAR PLUS HALF
    lamda1_ubarplus = ubarplus(i);
    lamda2_ubarplus = ubarplus(i)+a;
    lamda3_ubarplus = ubarplus(i)-a;
    
    lamda1plus_ubarplus = 0.5*(lamda1_ubarplus+abs(lamda1_ubarplus));
    lamda2plus_ubarplus = 0.5*(lamda2_ubarplus+abs(lamda2_ubarplus));
    lamda3plus_ubarplus = 0.5*(lamda3_ubarplus+abs(lamda3_ubarplus));

    lamda1minus_ubarplus = 0.5*(lamda1_ubarplus-abs(lamda1_ubarplus));
    lamda2minus_ubarplus = 0.5*(lamda2_ubarplus-abs(lamda2_ubarplus));
    lamda3minus_ubarplus = 0.5*(lamda3_ubarplus-abs(lamda3_ubarplus));
    %PUT THIS INTO BIG LAMDA PLUS AND BIG LAMDA MINUS FOR U PLUS HALF
    
    lamda_plus_u_plus(1,1) = lamda1plus_ubarplus;
    lamda_plus_u_plus(2,2) = lamda2plus_ubarplus;
    lamda_plus_u_plus(3,3) = lamda3plus_ubarplus;
    
    lamda_minus_u_plus(1,1) = lamda1minus_ubarplus;
    lamda_minus_u_plus(2,2) = lamda2minus_ubarplus;
    lamda_minus_u_plus(3,3) = lamda3minus_ubarplus;  
    
    
    
    
    %FOR UBAR MINUS HALF
    lamda1_ubarminus = ubarminus(i);
    lamda2_ubarminus = ubarminus(i)+a;
    lamda3_ubarminus = ubarminus(i)-a;   
    
    lamda1plus_ubarminus = 0.5*(lamda1_ubarminus+abs(lamda1_ubarminus));
    lamda2plus_ubarminus = 0.5*(lamda2_ubarminus+abs(lamda2_ubarminus));
    lamda3plus_ubarminus = 0.5*(lamda3_ubarminus+abs(lamda3_ubarminus));

    lamda1minus_ubarminus = 0.5*(lamda1_ubarminus-abs(lamda1_ubarminus));
    lamda2minus_ubarminus = 0.5*(lamda2_ubarminus-abs(lamda2_ubarminus));
    lamda3minus_ubarminus = 0.5*(lamda3_ubarminus-abs(lamda3_ubarminus));
    %PUT THIS INTO BIG LAMDA PLUS AND BIG LAMDA MINUS FOR U MINUS HALF
    
    lamda_plus_u_minus(1,1) = lamda1plus_ubarminus;
    lamda_plus_u_minus(2,2) = lamda2plus_ubarminus;
    lamda_plus_u_minus(3,3) = lamda3plus_ubarminus;
    
    lamda_minus_u_minus(1,1) = lamda1minus_ubarminus;
    lamda_minus_u_minus(2,2) = lamda2minus_ubarminus;
    lamda_minus_u_minus(3,3) = lamda3minus_ubarminus;
    
%% FIND X FOR U PLUS HALF 

    
    
     
    %HARD CODE THE X MATRIX    


    alpha = rho(i)./(a*sqrt(2.));
    x_matrix_uplus(1,1) = 1;
    x_matrix_uplus(1,2) = alpha;
    x_matrix_uplus(1,3) = alpha;
    x_matrix_uplus(2,1) = ubarplus(i);
    x_matrix_uplus(2,2) = alpha.*(ubarplus(i)+a);
    x_matrix_uplus(2,3) = alpha.*(ubarplus(i)-a);
    x_matrix_uplus(3,1) = 0.5*ubarplus(i).^2.;
    x_matrix_uplus(3,2) = alpha.*((0.5.*ubarplus(i).^2.)+(ubarplus(i).*a)+(a.^2./(gamma-1)));
    x_matrix_uplus(3,3) = alpha.*((0.5.*ubarplus(i).^2.)-(ubarplus(i).*a)+(a.^2./(gamma-1)));

     %FIND X INVERSE
     x_matrix_inverse_uplus = inv(x_matrix_uplus);

%% FIND X FOR U MINUS HALF 
     
     %HARD CODE THE X MATRIX
     
    alpha = rho(i)./(a*sqrt(2.));
    x_matrix_uminus(1,1) = 1;
    x_matrix_uminus(1,2) = alpha;
    x_matrix_uminus(1,3) = alpha;
    x_matrix_uminus(2,1) = ubarminus(i);
    x_matrix_uminus(2,2) = alpha.*(ubarminus(i)+a);
    x_matrix_uminus(2,3) = alpha.*(ubarminus(i)-a);
    x_matrix_uminus(3,1) = 0.5*ubarminus(i).^2;
    x_matrix_uminus(3,2) = alpha.*((0.5.*ubarminus(i).^2)+(ubarminus(i).*a)+(a.^2/(gamma-1)));
    x_matrix_uminus(3,3) = alpha.*((0.5.*ubarminus(i).^2)-(ubarminus(i).*a)+(a.^2/(gamma-1)));

     %FIND X INVERSE
     x_matrix_inverse_uminus = inv(x_matrix_uminus);
     
%% FIND ABARS

%FOR U PLUS HALF 
abarplus_uplus = x_matrix_uplus*lamda_plus_u_plus*x_matrix_inverse_uplus;
abarminus_uplus = x_matrix_uplus*lamda_minus_u_plus*x_matrix_inverse_uplus;

%FOR U MINUS HALF 

abarplus_uminus = x_matrix_uminus*lamda_plus_u_minus*x_matrix_inverse_uminus;
abarminus_uminus = x_matrix_inverse_uminus*lamda_minus_u_minus*x_matrix_inverse_uminus;
     
     
     
%% FIND FPLUS AND FMINUS


% %if split matters
fplus = abarplus_uplus*U_state(:,i)+ abarminus_uplus*U_state(:,i+1);
fminus = abarplus_uminus*U_state(:,i-1)+abarminus_uminus*U_state(:,i);

     
%if split doesnt matter

% fplus = abarplus_uplus*U_state(:,i)+ abarminus_uplus*U_state(:,i+1);
% fminus = abarplus_uplus*U_state(:,i-1)+abarminus_uplus*U_state(:,i);        


%% FIND U STATE AND UPDATE


%   U_state_update(:,i) = U_state(:,i) - (dt/dx).*(fplus-fminus);
%   U_state(:,i) = U_state_update(:,i);
%   
% 
%   t = t+ dt;
  end
  
  
  
 
  
  
  
  

end

 




