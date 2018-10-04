%============================%
%Shock tube problem(BTP) using FVS method by Steger_Warming splitting

clear all;
g=0;
gamma  = 1.4;
x_min=0;
x_max=2;
N=41;
d_x=(x_max-x_min)/N;
for i =  1:N                  % Initialization at t(time)=0
    if i<=N/2
        rho(i) = 2.0;    %density
        P(i) = 2.0;      %pressure
        u(i) = 0.0;      %velocity
    else
        rho(i) = 1.0;
        P(i) = 1.0;
        u(i) = 0.0;
       
    end
end
 E = P/(gamma-1)+0.5*rho.*u.^2; %energy
 U2 = rho.*u;       
 U3 = rho.*E;
t=0;
t_end = 12;
dt  = 0.01; 
while t<t_end
        
        a = (gamma*P./rho).^0.5;   % speed of sound
     
        L1 = u - a;                 %first eigen value
        L2 = u;                     %second eigen value
        L3 = u + a;                 %third eigen value
        
        L1P = 0.5*(L1+abs(L1)); 
        L2P = 0.5*(L2+abs(L2));
        L3P = 0.5*(L3+abs(L3));
        L1N = 0.5*(L1-abs(L1));
        L2N = 0.5*(L2-abs(L2));
        L3N = 0.5*(L3-abs(L3));
        
        H = .5*u.^2+a.^2/(gamma-1);         % H = 
%         (energy+pressure)/density
       
        FP1 = rho*0.5/gamma.*(L1P+2*(gamma-1)*L2P+L3P);  %first element of flux matrix (positive)
        FP2 = rho.*0.5/gamma.*((u-a).*L1P+2*(gamma-1)*u.*L2P+(u+a).*L3P);    %Second element of flux matrix
        FP3 = rho.*0.5/gamma.*((H-u.*a).*L1P+(gamma-1)*u.^2.*L2P+(H+u.*a).*L3P); %third element of flux matrix
        
        FN1 = rho.*0.5/gamma.*(L1N+2*(gamma-1)*L2N+L3N);
        FN2 = rho.*0.5/gamma.*((u-a).*L1N+2*(gamma-1)*u.*L2N+(u+a).*L3N);
        FN3 = rho.*0.5/gamma.*((H-u.*a).*L1N+(gamma-1)*u.^2.*L2N+(H+u.*a).*L3N);
        
        for i = 1:N-1                     %intercell numerical flux
            Fhp1(i) = FP1(i)+FN1(i+1);
            Fhn1(i+1) = FP1(i)+FN1(i+1);
            Fhp2(i) = FP2(i)+FN2(i+1);
            Fhn2(i+1) = FP2(i)+FN2(i+1);
            Fhp3(i) = FP3(i)+FN3(i+1);
            Fhn3(i+1) = FP3(i)+FN3(i+1);
        end
  
       for i=2:N-1
        rhon(i) = rho(i)-dt*(Fhp1(i)-Fhn1(i));        % Density at t = t+dt
        U2(i) = rho(i).*u(i)-dt*(Fhp2(i)-Fhn2(i));    % U2 at t=t+dt       
        U3(i) = rho(i).*E(i)-dt*(Fhp3(i)-Fhn3(i));    % U3 at t = t+dt
       end
       
      % Boundary Conditions
       rhon(N) = rhon(N-1);
       U2(N) = U2(N-1);
       U3(N) = U3(N-1);
       rhon(1) = rhon(2);
       U2(1) = U2(2);
       U3(1) = U3(2);
       
       u = U2./rhon;                       % velocity at t+dt
       E  = U3./rhon;                      % energy at t+dt
       P = (gamma-1)*(E-0.5*rhon.*u.^2);   %pressure at t+dt
       
       rho = rhon; % new density       
       t = t+dt;          % time increment                    
end
% plot(1:N,P,'LineWidth',2)
%         title('Pressure vs position');
%         xlabel('position');
%         ylabel('Pressure');
%         grid on;
%         figure
% plot(1:N,rho,'LineWidth',2)
%         title('Density vs position');
%         xlabel('position');
%         ylabel('Density');
       
plot(1:N,u,'LineWidth',2)
        title('Velocity vs position');
        xlabel('position');
        ylabel('Velocity');