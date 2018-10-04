function [F] = build_Flux( U )
gamma=1.4;
u1=U(1,:);
u2=U(2,:);
u3=U(3,:);
F(1,:)=u2;
F(2,:)=u2.^2/u1+(gamma-1).*(u3-u2.^2./(2.*u1));
F(3,:)=u2./u1.*(u3+(gamma-1).*(u3-u2.^2./(2.*u1)));
end

