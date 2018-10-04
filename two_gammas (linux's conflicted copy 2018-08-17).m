clear all
clc
% clf
%% INITIALIZATION 
%BASIC PARAMETERS
p4 = 2.0; 
p1 = 1.0; 
rho4 = 2.0; 
rho1 = 1.0; 
imax = 9; 
xmin = 0; 
xmax = 2.; 
dx = 2/(imax-1); 

global rho_vector
global p_vector
global velocity_vector
global dtnew


%SET UP VECTORS

x = 0:dx:2;
x0 = find(x==1.0);

u = zeros(1,imax);

%SET UP STATE VECTORS
USTATE = zeros(3,imax);
USTATE_UPDATE = zeros(3,imax);

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




timestep = 1;


% maxtimestep = 2;

% maxtimestep = input('Enter max time step: ');



%% DECLARE PHI 



phi = zeros(1,imax);
phi2 = 1.0;
phi1 = 0.0;

for i = 1:imax
   if (x(i)>=1)
       phi(i) = phi1;
   elseif (x(i)<1)
       phi(i) = phi2;
   end
    
end




phi_next = zeros(1,imax);

phi_initial  = phi;








%% CALCULATING GAMMA USING Chaklas Eq 19b 

gamma2 = 1.6;
gamma1 = 1.4;


gamma_at_i = zeros(1,imax);

for i = 1:imax
    temp = (phi(i)./(gamma1-1)) + ((1-phi(i))./(gamma2-1));
    gamma_at_i(i) = 1 + 1./temp;
end

               


%% Chaklas Eq 19c

%assume P1_infty = p1, and P2_infty = p4



t_ambient = 300;
R = 287;
p1_infty = rho1*R*t_ambient;
p2_infty = rho4*R*t_ambient;

p_infty = ((gamma_at_i-1)./(gamma_at_i)).*(((phi.*gamma1*p1_infty)./(gamma1-1))  + ...
(((1-phi).*(gamma2*p2_infty))/(gamma2-1)));

gamma_initial = gamma_at_i;






%% INITIAL CONDITION
% RHO
USTATE(1,1:x0) = rho4;
USTATE(1,x0:imax) = rho1;

% FOR E 
USTATE(3,1:x0) = p4./(gamma_at_i(1:x0)-1);
USTATE(3,x0:imax) = p1./(gamma_at_i(x0:imax)-1);

USTATE_0 = USTATE;
USTATE_UPDATE = USTATE; 



% %PRESSURE 
% 
% PRESSSURE_INITIAL(1:x0) = (gamma_at_i(1:x0)-1).*USTATE(3,1:x0);
% PRESSSURE_INITIAL(x0:imax) = (gamma_at_i(x0:imax)-1).*USTATE(3,x0:imax);



%PRESURE (using Chaklas Eq 19a)

PRESSSURE_INITIAL(1:x0) = ((gamma_at_i(1:x0)-1).*USTATE(3,1:x0))  - ...
(gamma_at_i(1:x0).*p_infty(1:x0));
PRESSSURE_INITIAL(x0:imax) = ((gamma_at_i(x0:imax)-1).*USTATE(3,x0:imax)) - ...
(gamma_at_i(x0:imax).*p_infty(x0:imax));


PRESSURE = zeros(1,imax);






% PRESSURE = PRESSSURE_INITIAL;

%% MAIN LOOP
for timestep = 1: 1
    
    
    
    %% PART 1: CALCULATE DT 
 
  
  for i = 1:imax
     u_at_i = USTATE_UPDATE(2,:)./USTATE_UPDATE(1,:); %%rho u by rho 
     rho_at_i = USTATE_UPDATE(1,:);
     e_at_i = USTATE_UPDATE(3,:);
     p_at_i = (e_at_i-0.5.*rho_at_i.*u_at_i.^2);
     a_at_i = sqrt(gamma_at_i(i)*p_at_i./rho_at_i);
     abs_u_plus_a_at_i= abs(u_at_i+a_at_i);
     dt_at_i = dx./abs_u_plus_a_at_i;%this is an array of dt at i
     dt_smallest = min(dt_at_i); %smallest dt in dt array at i
     real_dt = 0.9*dt_smallest;
     
    
    
  end
  
  dt = real_dt;
  
   %% PART 2: DEAL WITH U PLUS HALF
   for i = 2:imax-1
       
 
       
    USTATE_PLUS(:,:) = 0.5*(USTATE(:,i)+USTATE(:,i+1));
    
    
    rho_plus = USTATE_PLUS(1,:);   
   
    u_plus = USTATE_PLUS(2,:)./USTATE_PLUS(1,:);
    
    
    
    p_plus = (USTATE_PLUS(3,:).*(gamma_at_i(i)-1)) - (0.5).*(gamma_at_i(i)-1).*rho_plus.*u_plus^2;   
    
  
    a_plus = sqrt(gamma_at_i(i).*p_plus./rho_plus);
    
    
    
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
    Ca_i_plus(1,1) = 1.0;
    Ca_i_plus(1,2) = 0.0;
    
    Ca_i_plus(1,3) = -1./(a_plus.^2);
    Ca_i_plus(2,1) = 0.0;
    Ca_i_plus(2,2) = rho_plus.*a_plus;
    Ca_i_plus(2,3) = 1.0;
    
    Ca_i_plus(3,1) = 0.0;
    Ca_i_plus(3,2) = -rho_plus.*a_plus;
    Ca_i_plus(3,3) = 1.0;
    
    
      
    %Form S plus     
    beta = gamma_at_i(i)-1;
    alpha_plus = (u_plus.^2)./2;
    S_i_plus(1,1) = 1.0;
    S_i_plus(1,2) = 0.0;
    S_i_plus(1,3) = 0.0;
    S_i_plus(2,1) = -u_plus./rho_plus;
    S_i_plus(2,2) = 1.0./rho_plus;
    S_i_plus(2,3) = 0.0;
    S_i_plus(3,1) = alpha_plus.*beta;
    S_i_plus(3,2) = -u_plus.*beta;
    S_i_plus(3,3) = beta;
    
    
    %Form Ca inverse plus 
    
    Ca_inverse_plus(1,1) = 1.0;
    Ca_inverse_plus(1,2) = 1.0./(2.*a_plus.^2);
    Ca_inverse_plus(1,3) = 1.0./(2.*a_plus.^2);
    
    Ca_inverse_plus(2,1) = 0.0;
    
    Ca_inverse_plus(2,2) = 1.0./(2.*rho_plus.*a_plus);
    Ca_inverse_plus(2,3) = -1.0./(2.*rho_plus.*a_plus);
    
    Ca_inverse_plus(3,1) = 0.0;
    Ca_inverse_plus(3,2) = 0.5;
    Ca_inverse_plus(3,3) = 0.5;


    
    %Form S inverse plus 
    
    
    S_inverse_plus(1,1) = 1.0;
    S_inverse_plus(1,2) = 0.0;
    S_inverse_plus(1,3) = 0.0;
    S_inverse_plus(2,1) = u_plus;
    S_inverse_plus(2,2) = rho_plus;
    S_inverse_plus(2,3) = 0.0;
    S_inverse_plus(3,1) = alpha_plus;
    S_inverse_plus(3,2) = rho_plus.*u_plus;
    S_inverse_plus(3,3) = 1.0./beta;
    
   
  
   
    Abarplus_plus(:,:) = S_inverse_plus(:,:)*Ca_inverse_plus(:,:)*lamda_plus_i_plus(:,:)*Ca_i_plus(:,:)*S_i_plus(:,:);

% %     
%     Abarplus_plus(:,:) = S_inverse_plus(:,:)*Ca_inverse_plus(:,:);
%     Abarplus_plus(:,:) = Abarplus_plus(:,:)*lamda_plus_i_plus(:,:);
%     Abarplus_plus(:,:) = Abarplus_plus(:,:)*Ca_i_plus(:,:);
%     Abarplus_plus(:,:) = Abarplus_plus(:,:)*S_i_plus(:,:);

    Abarminus_plus(:,:) = S_inverse_plus(:,:)*Ca_inverse_plus(:,:)*lamda_minus_i_plus(:,:)*Ca_i_plus(:,:)*S_i_plus(:,:);

    
    FPLUS(:,i) = Abarplus_plus(:,:)*USTATE(:,i)+Abarminus_plus(:,:)*USTATE(:,i+1);

   end
    
    
    
    
    
    %% PART 3: DEAL WITH U MINUS HALF
    
    for i = 2:imax-1
    
    USTATE_MINUS(:,:) = 0.5*(USTATE(:,i)+USTATE(:,i-1));
    
    %PULL RHO, P, AT MINUS HALF
    
    rho_minus = USTATE_MINUS(1,:);
   
  
    
    u_minus = USTATE_MINUS(2,:)./USTATE_MINUS(1,:);
    p_minus = (USTATE_MINUS(3,:).*(gamma_at_i(i)-1)) - (0.5)*(gamma_at_i(i)-1)*rho_minus.*u_minus.^2;    

    a_minus = sqrt(gamma_at_i(i)*p_minus/rho_minus);
    
    
    
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
    
    
    
    %FORMING Ca, S for I MINUS HALF 
    %Form Ca minus
    Ca_i_minus(1,1) = 1.0;
    Ca_i_minus(1,3) = -1./(a_minus.^2);
    Ca_i_minus(2,2) = rho_minus.*a_minus;
    Ca_i_minus(2,3) = 1.0;
    Ca_i_minus(3,2) = -rho_minus.*a_minus;
    Ca_i_minus(3,3) = 1.0;
    %Form S minus     
    beta = gamma_at_i(i)-1;
    alpha_minus = (u_minus.^2)./2;
    S_i_minus(1,1) = 1.0;
    S_i_minus(2,1) = -u_minus./rho_minus;
    S_i_minus(2,2) = 1.0./rho_minus;
    S_i_minus(3,1) = alpha_minus.*beta;
    S_i_minus(3,2) = -u_minus.*beta;
    S_i_minus(3,3) = beta;
    
    %Form Ca inverse minus 
    
    Ca_inverse_minus(1,1) = 1.0;
    Ca_inverse_minus(1,2) = 1.0./(2.*a_minus.^2);
    Ca_inverse_minus(1,3) = 1.0./(2.*a_minus.^2);
    Ca_inverse_minus(2,2) = 1.0./(2.*rho_minus.*a_minus);
    Ca_inverse_minus(2,3) = -1.0./(2.*rho_minus.*a_minus);
    Ca_inverse_minus(3,2) = 0.5;
    Ca_inverse_minus(3,3) = 0.5;
    
    %Form S inverse minus 
    
    
    S_inverse_minus(1,1) = 1.0;
    S_inverse_minus(2,1) = u_minus;
    S_inverse_minus(2,2) = rho_minus;
    S_inverse_minus(3,1) = alpha_minus;
    S_inverse_minus(3,2) = rho_minus.*u_minus;
    S_inverse_minus(3,3) = 1.0./beta;    
    
    
    
    Abarplus_minus(:,:) = S_inverse_minus(:,:)*Ca_inverse_minus(:,:)*lamda_plus_i_minus(:,:)*Ca_i_minus(:,:)*S_i_minus(:,:);

    Abarminus_minus(:,:) = S_inverse_minus(:,:)*Ca_inverse_minus(:,:)*lamda_minus_i_minus(:,:)*Ca_i_minus(:,:)*S_i_minus(:,:);
    
    FMINUS(:,i) = Abarplus_minus(:,:)*USTATE(:,i-1)+Abarminus_minus(:,:)*USTATE(:,i);
    
    end
    
   
    
    %% PART 4: FINITE DIFFERENCE EQUATION
    
    for i = 2:imax-1    
%         USTATE_UPDATE(:,i) = USTATE(:,i) - (dt/dx)*(FPLUS(:,i)-FMINUS(:,i)); 
        USTATE_UPDATE(:,i) = (dt/dx)*(FPLUS(:,i)-FMINUS(:,i)); 
             
    end
   
  
  
  %% SET BC AND PLOT VARIABLES FOR NUMERICAL 

  
  USTATE_UPDATE(:,imax) = USTATE_UPDATE(:,imax-1);
  USTATE = USTATE_UPDATE; 
  RHO_E = USTATE(3,:); 
  RHO = USTATE(1,:);
  
  VELOCITY = USTATE(2,:)./RHO;

%   %CALCULATE PRESSURE, not using 19a, 19c  
  for i = 1:imax
    PRESSURE(i) = (gamma_at_i(i)-1)*(RHO_E(i)-0.5.*RHO(i)*VELOCITY(i)^2);
  end
%  
%     



%  %% Chaklas 19a and c 
%  
%  
% p_infty = ((gamma_at_i-1)./(gamma_at_i)).*(((phi.*gamma1*p1_infty)./(gamma1-1))  + ...
% (((1-phi).*(gamma2*p2_infty))/(gamma2-1)));
% 
% PRESSURE = ((gamma_at_i-1).*USTATE(3,:)) - ...
% (gamma_at_i.*p_infty);

 
 



    
  %% PART 5: ANALYTICAL PART  (call analytical, use updated dt);
  
  for i = 1:imax
    max_ANALYTICAL_shock_tube(dt,gamma_at_i(i));
  end
  
  
  
  %% PART 6: CALCULATE PHI
     
    phi_next = phi;
    
    %Loop and find Uminus, Uplus in the VELOCITY
    
    for i = 2:imax-2
        uplusDCU = max(0,VELOCITY(1,:));
        uminusDCU = min(0,VELOCITY(1,:));
    
        %the finite difference expressions
        forward_diffX_PHI = phi(i+1)-phi(i);
        backward_diffX_PHI = phi(i)-phi(i-1);

        
        %main equation            
        phi_next(i)= phi(i)-(dt/dx)*(uplusDCU(i)*backward_diffX_PHI+uminusDCU(i)*...
        forward_diffX_PHI); 
               
   
    end
    
    
    phi = phi_next;
 
      

    
    %% PART 7: RE CALCULATE GAMMA AT I 
    
%     rhs_19b_updated = (phi./(gamma1-1)) + ((1-phi)./(gamma2-1));
%     gamma_at_i = 1 + 1./(rhs_19b_updated);
% 
       
    for i = 1:imax
        temp = (phi(i)./(gamma1-1)) + ((1-phi(i))./(gamma2-1));
        gamma_at_i(i) = 1 + 1./temp;
    end


                        



    
    
    %% PART 8: PLOTTING ANIMATION
   
    
    
     disp(timestep)
     figure(1)
     clf()
     subplot(2,3,1)
     plot(x,PRESSURE,'-xr');
     hold on
     grid on
     plot(x,p_vector,'-ob');
     legend({'Steger Warming','Analytical'},'FontSize',12,'Location','northeast');
     xlabel('X','FontSize',18);
     title('X vs PRESSURE STEGER-WARMING EXPLICIT','FontSize',10);
     ylabel('PRESSURE','FontSize',12); 
     xt = get(gca, 'XTick');
     set(gca, 'FontSize', 12)
 
     subplot(2,3,2) 
     plot(x,RHO,'-xr');
     ylim([0.5 2]);
     hold on
     grid on
     plot(x,rho_vector,'-ob');
     legend({'Steger Warming','Analytical'},'FontSize',12,'Location','southwest');
     xlabel('X','FontSize',12);
     title('X vs DENSITY STEGER-WARMING EXPLICIT','FontSize',12)
     ylabel('DENSITY','FontSize',12);
     xt = get(gca, 'XTick');
     set(gca, 'FontSize', 12);
     
         
     subplot(2,3,3) 
     plot(x,VELOCITY,'-xr');
     hold on
     grid on
     plot(x,velocity_vector,'-ob');
     legend({'Steger Warming','Analytical'},'FontSize',12,'Location','northwest');
     xlabel('X','FontSize',12);
     title('X vs VELOCITY STEGER-WARMING EXPLICIT','FontSize',12)
     ylabel('VELOCITY','FontSize',12);
     xt = get(gca, 'XTick');
     set(gca, 'FontSize', 12);
     
     
     
 
     subplot(2,3,4)
     plot(x,phi_initial,'-ob')
     hold on
     plot(x,phi_next,'-xr')
     title('Species Phi Interface','FontSize',12)
     xlabel('X','FontSize',12);
     ylabel('PHI','FontSize',12);
     xt = get(gca, 'XTick');
     set(gca, 'FontSize', 12)
     legend('Initial','After')
 
 
     
     
     subplot(2,3,5)
     plot(x,gamma_initial,'-ob')
     hold on
     plot(x,gamma_at_i,'-xr')
     title("Gamma Profile", "FontSize",12);
     xlabel("X","FontSize",12)
     ylabel('GAMMA','FontSize',12);
     xt = get(gca, 'XTick');
     set(gca, 'FontSize', 12)
     legend('Initial','After')
     pause(0.01)
%     
    
    
    
    timestep = timestep + dt;
    
    
    
end 














