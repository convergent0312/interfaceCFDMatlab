%------------------------------------------------------------------------%
%---------Sod Shock tube by Richtmyer method-----------------------------%
%------------------------------------------------------------------------%


clc
clear all
close all

%---initial parameters----------------------------------------------------

nx = 81;
dx = .25;
dt = .0002  ; 
nt=50;
x=-10:0.25:10;
U(1:3,1:81)=1;

%----initial condition of left side of the tube----------------------------

rho=1;u=0;p=100000;
gamma = 1.6;
e = p / ((gamma-1)*rho);
eT = e + u^2/2;
u=[rho,rho*u,rho*eT];

for i=1:3
U(i,1:40)=u(i);
end

%----initial condition of right side of the tube-----------------------------


rho=0.125;u=0;p=10000;
gamma = 1.4;
e = p / ((gamma-1)*rho);
eT = e + u^2/2;
u=[rho,rho*u,rho*eT];

for i=1:3
U(i,41:81)=u(i);
end

%-------Richtmyer method implementation------------------------------------

UN=ones(3,nx);
UN_plus=ones(3,nx);
UN_minus=ones(3,nx);

for i=1:nt
    UN_plus(:,1:80)=0.5*(U(:,2:81)+U(:,1:80))-dt/(2*dx)*...
        (build_Flux(U(:,2:81))-build_Flux(U(:,1:80)));
     UN_minus(:,2:81) = UN_plus(:,1:80);  
        UN(:,2:80) = U(:,2:80) - dt/dx *...
        (build_Flux(UN_plus(:,2:80)) - build_Flux(UN_minus(:,2:80)));
     UN(:,1)=UN(:,2);
     UN(:,81)=UN(:,80);
     U=UN;
end

%---Calculating velocity, pressure and density---------------------------

vel = U(2,:)./U(1,:);
pres = (gamma - 1)*(U(3,:) - .5 * U(2,:).^2 ./ U(1,:));
rho = U(1,:);

%---plotting the results-------------------------------------------------%

figure(1)
subplot(131)
plot(x,vel,'linewidth',2)
   h=gca; 
   get(h,'FontSize') 
   set(h,'FontSize',14)
   xlabel('X','fontSize',14);
   ylabel('Velocity','fontSize',14);
   subplot(132)
   plot(x,pres,'linewidth',2) 
   h=gca; 
   get(h,'FontSize') 
   set(h,'FontSize',14)
   xlabel('X','fontSize',14);
   ylabel('Pressure','fontSize',14);
   subplot(133)
   plot(x,rho,'linewidth',2) 
   h=gca; 
   get(h,'FontSize') 
   set(h,'FontSize',14)
   xlabel('X','fontSize',14);
   ylabel('density','fontSize',14);
   fh = figure(1);
   set(fh, 'color', 'white'); 
   
   %----------------------------------------------------------%