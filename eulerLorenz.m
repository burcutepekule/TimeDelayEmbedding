function [x,y,z,t] = eulerLorenz( sig, beta, rho, IC, numPoints, dt )
x(1)=IC(1);
y(1)=IC(2);
z(1)=IC(3);
t(1)=0;
for k=1:numPoints % Time loop
    fx=sig*(y(k)-x(k)); % RHS of x equation
    fy=-x(k)*z(k)+rho*x(k)-y(k); % RHS of y equation
    fz=x(k)*y(k)-beta*z(k); % RHS of z equation
    x(k+1)=x(k)+dt*fx; % Find new x
    y(k+1)=y(k)+dt*fy; % Find new y
    z(k+1)=z(k)+dt*fz; % Find new z
    t(k+1)=t(k)+dt; % Find new t
end % Close time loop
end

