clear all 
clc

ustate = zeros(3,9);

ustate(1,1:5) = 10;
ustate(1,5:9) = 4;

for n = 1:1
    disp(ustate(1,:));
    for i = 2:8
        umid = 0.5*(ustate(:,i)+ustate(:,i+1));
    end
    disp(ustate(1,:));
    
end
