clear 
clc
clf
%% INITIALIZATION 
%BASIC PARAMETERS
gamma = 1.4;
p4 = 2.0; 
p1 = 1.0; 
rho4 = 2.0; 
rho1 = 1.0; 
imax = 41; 
xmin = 0; 
xmax = 2.; 
dx = 2/(imax-1); 

global rho_vector
global p_vector
global velocity_vector


%SET UP VECTORS

x = 0:dx:2;
x0 = find(x==1.0);

u = zeros(1,imax);

%SET UP STATE VECTORS
USTATE = zeros(3,imax);
% USTATE_UPDATE = zeros(3,imax);

%SET UP MATRICES 

lamda_plus_i_plus = zeros(3,3);
lamda_minus_i_plus = zeros(3,3);
Ca_i_plus = zeros(3,3);
Ca_inverse_plus = zeros(3,3);
S_i_plus = zeros(3,3);
S_inverse_plus = zeros(3,3);
FPLUS = zeros(3,imax);
xplus = zeros(3,3);
xplusinv = zeros(3,3);
xminus= zeros(3,3);
xminusinv= zeros(3,3);
% FPLUS = zeros(3,3);


lamda_plus_i_minus = zeros(3,3);
lamda_minus_i_minus = zeros(3,3);
Ca_i_minus = zeros(3,3);
Ca_inverse_minus = zeros(3,3);
S_i_minus = zeros(3,3);
S_inverse_minus = zeros(3,3);
FMINUS = zeros(3,imax);
% FMINUS = zeros(3,3);


Abarplus_plus= zeros(3,3);
Abarminus_plus= zeros(3,3);

Abarplus_minus= zeros(3,3);
Abarminus_minus= zeros(3,3);


% INITIAL CONDITIONS
% RHO
USTATE(1,1:x0) = rho4;
USTATE(1,x0:imax) = rho1;

% FOR E 
USTATE(3,1:x0) = p4/(gamma-1);
USTATE(3,x0:imax) = p1/(gamma-1);

USTATE_UPDATE = USTATE; 

timestep = 0;
maxtimestep = 18;





%% MAIN LOOP
while timestep <maxtimestep
  %% PART 1: CALCULATE DT 
  
  for i = 2:imax-1
     u_at_i = USTATE_UPDATE(2,:)./USTATE_UPDATE(1,:); %%rho u by rho 
     rho_at_i = USTATE_UPDATE(1,:);
     e_at_i = USTATE_UPDATE(3,:);
     p_at_i = (e_at_i-0.5.*rho_at_i.*u_at_i.^2);
     a_at_i = sqrt(gamma*p_at_i./rho_at_i);
     abs_u_plus_a_at_i= abs(u_at_i+a_at_i);
     dt_at_i = dx./abs_u_plus_a_at_i;
     dt_smallest = min(dt_at_i); %smallest dt in dt at i
     real_dt = 0.9*dt_smallest;
     
    
    
  end
  
  dt = real_dt;
  dt = 0.02;
  
   %% PART 2: DEAL WITH U PLUS HALF
    for i = 2:imax-1
    USTATE_PLUS(:,:) = 0.5*(USTATE(:,i)+USTATE(:,i+1));
    
    
    rho_plus = USTATE_PLUS(1,:);   
   
    u_plus = USTATE_PLUS(2,:)./USTATE_PLUS(1,:);
    
    p_plus = (USTATE_PLUS(3,:).*(gamma-1)) - (0.5).*(gamma-1).*rho_plus.*u_plus^2;    
  
    a_plus = sqrt(gamma.*p_plus./rho_plus);
    
    %FORM LAMBDA PLUS AND MINUS FOR I PLUS HALF 
    lamda1plus = u_plus;
    lamda2plus = u_plus+a_plus;
    lamda3plus = u_plus-a_plus;
    
    lamda1plus_plus = 0.5*(lamda1plus+abs(lamda1plus));
    lamda2plus_plus = 0.5*(lamda2plus+abs(lamda2plus));
    lamda3plus_plus = 0.5*(lamda3plus+abs(lamda3plus));
    
    lamda1minus_plus = 0.5*(lamda1plus-abs(lamda1plus));
    lamda2minus_plus = 0.5*(lamda2plus-abs(lamda2plus));
    lamda3minus_plus = 0.5*(lamda3plus-abs(lamda3plus));
    
    %LAMDA PLUS FOR I PLUS HALF 
    lamda_plus_i_plus(1,1) = lamda1plus_plus;
    lamda_plus_i_plus(2,2) = lamda2plus_plus;
    lamda_plus_i_plus(3,3) = lamda3plus_plus;

    
    %LAMDA PLUS FOR I MINUS HALF 
    lamda_minus_i_plus(1,1) = lamda1minus_plus;
    lamda_minus_i_plus(2,2) = lamda2minus_plus;
    lamda_minus_i_plus(3,3) = lamda3minus_plus;
    
    
    
    %FORMING Ca, S for I PLUS HALF 
    %Form Ca plus
    Ca_i_plus(:,:) = [1.0 0.0 -1./(a_plus.^2); 0.0 rho_plus.*a_plus 1.0; 0.0 -rho_plus.*a_plus 1.0]; 
    
    
    
      
    %Form S plus     
    beta = gamma-1;
    alpha_plus = (u_plus.^2)./2;
    S_i_plus(:,:) = [1.0 0.0 0.0; -u_plus./rho_plus 1.0./rho_plus 0.0; alpha_plus.*beta -u_plus.*beta beta];
    
    
    %Form Ca inverse plus 
    
    Ca_inverse_plus(:,:) = [1.0 1.0./(2.*a_plus.^2) 1.0./(2.*a_plus.^2); 0.0 1.0./(2.*rho_plus.*a_plus) -1.0./(2.*rho_plus.*a_plus); 0.0 0.5 0.5];
    


    
    %Form S inverse plus 
    
    
    S_inverse_plus(:,:) = [1.0 0.0 0.0; u_plus rho_plus 0.0; alpha_plus rho_plus.*u_plus 1.0./beta];
    
    
    
  

    
    
  
   
    Abarplus_plus(:,:) = S_inverse_plus(:,:)*Ca_inverse_plus(:,:)*lamda_plus_i_plus(:,:)*Ca_i_plus(:,:)*S_i_plus(:,:);

    Abarminus_plus(:,:) = S_inverse_plus(:,:)*Ca_inverse_plus(:,:)*lamda_minus_i_plus(:,:)*Ca_i_plus(:,:)*S_i_plus(:,:);

    
    FPLUS(:,i) = Abarplus_plus(:,:)*USTATE(:,i)+Abarminus_plus(:,:)*USTATE(:,i+1);

    
    
    
    
    %% PART 3: DEAL WITH U MINUS HALF
    
    
    USTATE_MINUS(:,:) = 0.5*(USTATE(:,i)+USTATE(:,i-1));
    
    %PULL RHO, P, AT MINUS HALF
    
    rho_minus = USTATE_MINUS(1,:);
    u_minus = USTATE_MINUS(2,:)./USTATE_MINUS(1,:);
    p_minus = (USTATE_MINUS(3,:).*(gamma-1)) - (0.5)*(gamma-1)*rho_minus.*u_minus.^2;    

    a_minus = sqrt(gamma*p_minus/rho_minus);
    
    %FORM LAMBDA PLUS AND MINUS FOR I MINUS HALF 
    lamda1minus = u_minus;
    lamda2minus = u_minus+a_minus;
    lamda3minus = u_minus-a_minus;
    
    lamda1plus_minus = 0.5*(lamda1minus+abs(lamda1minus));
    lamda2plus_minus = 0.5*(lamda2minus+abs(lamda2minus));
    lamda3plus_minus = 0.5*(lamda3minus+abs(lamda3minus));
    
    lamda1minus_minus = 0.5*(lamda1minus-abs(lamda1minus));
    lamda2minus_minus = 0.5*(lamda2minus-abs(lamda2minus));
    lamda3minus_minus = 0.5*(lamda3minus-abs(lamda3minus));
    
    %LAMDA PLUS FOR I MINUS HALF 
    lamda_plus_i_minus(1,1) = lamda1plus_minus;
    lamda_plus_i_minus(2,2) = lamda2plus_minus;
    lamda_plus_i_minus(3,3) = lamda3plus_minus;
    
    %LAMDA MINUS FOR I MINUS HALF 
    lamda_minus_i_minus(1,1) = lamda1minus_minus;
    lamda_minus_i_minus(2,2) = lamda2minus_minus;
    lamda_minus_i_minus(3,3) = lamda3minus_minus;
    
    
    
       %FORMING Ca, S for I PLUS HALF 
    %Form Ca plus
    Ca_i_minus(:,:) = [1.0 0.0 -1./(a_minus.^2); 0.0 rho_minus.*a_minus 1.0; 0.0 -rho_minus.*a_minus 1.0]; 
    
    
    
      
    %Form S plus     
    beta = gamma-1;
    alpha_minus = (u_minus.^2)./2;
    S_i_minus(:,:) = [1.0 0.0 0.0; -u_minus./rho_minus 1.0./rho_minus 0.0; alpha_minus.*beta -u_minus.*beta beta];
    
    
    %Form Ca inverse plus 
    
    Ca_inverse_minus(:,:) = [1.0 1.0./(2.*a_minus.^2) 1.0./(2.*a_minus.^2); 0.0 1.0./(2.*rho_minus.*a_minus) -1.0./(2.*rho_minus.*a_minus); 0.0 0.5 0.5];
    


    
    %Form S inverse plus 
    
    
    S_inverse_minus(:,:) = [1.0 0.0 0.0; u_minus rho_minus 0.0; alpha_minus rho_minus.*u_minus 1.0./beta];
    
    
    
    Abarplus_minus(:,:) = S_inverse_minus(:,:)*Ca_inverse_minus(:,:)*lamda_plus_i_minus(:,:)*Ca_i_minus(:,:)*S_i_minus(:,:);

    Abarminus_minus(:,:) = S_inverse_minus(:,:)*Ca_inverse_minus(:,:)*lamda_minus_i_minus(:,:)*Ca_i_minus(:,:)*S_i_minus(:,:);
    
    FMINUS(:,i) = Abarplus_minus(:,:)*USTATE(:,i-1)+Abarminus_minus(:,:)*USTATE(:,i);
    
    %% PART 4: FINITE DIFFERENCE EQUATION
    
    USTATE_UPDATE(:,i) = USTATE(:,i) - (dt/dx)*(FPLUS(:,i)-FMINUS(:,i)); 
%     USTATE(:,i) = USTATE_UPDATE(:,i);

     
    end
   
  
  
  %% SET BC AND PLOT VARIABLES FOR NUMERICAL 
  
  USTATE_UPDATE(:,imax) = USTATE_UPDATE(:,imax-1);
  USTATE = USTATE_UPDATE;
  
%   USTATE_UPDATE(1,imax) = USTATE_UPDATE(1,imax-1);
%   USTATE_UPDATE(2,imax) = USTATE_UPDATE(2,imax-1);
%   USTATE_UPDATE(3,imax) = USTATE_UPDATE(3,imax-1);
  
%   USTATE_UPDATE(1,1) = USTATE_UPDATE(1,2);
%   USTATE_UPDATE(2,1) = USTATE_UPDATE(2,2);
%   USTATE_UPDATE(3,1) = USTATE_UPDATE(3,2);
   
  ENERGY = USTATE_UPDATE(3,:); 
  RHO = USTATE_UPDATE(1,:);
  
  VELOCITY = USTATE_UPDATE(2,:)./RHO;
  PRESSURE = ((gamma-1)*ENERGY-(gamma-1)*0.5.*(((USTATE_UPDATE(2,:)).^2)./(RHO)));
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% PART 5: ANALYTICAL PART  (call analytical, use dt);
  
  shocktube(dt);
  timestep = timestep+1;

end 

  
  
  
  
  
  
    
  figure(1)
  xlabel('X');
  ylabel('PRESSURE');
  title('X vs PRESSURE')
  plot(x,PRESSURE);
  hold on
  plot(x,p_vector);
  legend('Steger Warming','Analytical');
%   
%   
%   xlabel('X');
%   ylabel('DENSITY');
%   title('X vs DENSITY')
%   plot(x,RHO);
%   hold on
%   plot(x,rho_vector);
%   legend('Steger Warming','Analytical');
  
%   
 
%   xlabel('X');
%   ylabel('VELOCITY');
%   title('X vs VELOCITY')
%   plot(x,VELOCITY);
%   hold on
%   plot(x,velocity_vector);
%   legend('Steger Warming','Analytical');
%   
  
  










