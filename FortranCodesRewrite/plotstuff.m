clear all 
clc

% LOADING DATA
load('rdata.dat')
load('prdata.dat')
load('vdata.dat')
load('phidata.dat')
load('gammadata.dat')



figure(1)
subplot(2,3,1)
plot(rdata(:,1),rdata(:,2), '-xr')
hold on
plot(rdata(:,1),rdata(:,3), '-ob')
title('RHO')
legend('INITIAL', 'AFTER')
xlabel('X')
ylabel('RHO')

subplot(2,3,2)

plot(prdata(:,1),prdata(:,2), '-xr')
hold on
plot(prdata(:,1),prdata(:,3), '-ob')
title('PRESSURE')
legend('INITIAL', 'AFTER')
xlabel('X')
ylabel('PRESSURE')

subplot(2,3,3)

plot(vdata(:,1),vdata(:,2), '-xr')
hold on
plot(vdata(:,1),vdata(:,3), '-ob')
title('VELOCITY')
legend('INITIAL', 'AFTER')
xlabel('X')
ylabel('VELOCITY')

subplot(2,3,4)

plot(phidata(:,1),phidata(:,2), '-xr')
hold on
plot(phidata(:,1),phidata(:,3), '-ob')
title('PHI')
legend('INITIAL', 'AFTER')
xlabel('X')
ylabel('PHI')


subplot(2,3,5)
plot(gammadata(:,1),gammadata(:,2), '-xr')
hold on
plot(gammadata(:,1),gammadata(:,3), '-ob')
title('GAMMA')

legend('INITIAL', 'AFTER')
xlabel('X')
ylabel('GAMMA')



