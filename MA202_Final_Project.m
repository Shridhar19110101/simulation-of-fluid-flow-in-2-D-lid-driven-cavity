%%%solving vorticity/stream function and velocity function
%%for lid driven cavity
clear; 
close all
%declaration of terms
N = 30;
l = 1; %length of grid
U_wall = 1; %Velocity of moving wall
rho = 1; mu = 0.01; % Density; Dynamic Viscosity;
dt = 0.001;% Time Step
maxIt = 50000; maxe = 1e-7; % Max iter; Max error
% SETUP 2D GRID
h=l/(N-1); %step size of grid 
i = 2:N-1; 
j = 2:N-1; 
% allocating MATRIXES
w = zeros(N,N);  
psi = zeros(N,N);
wp = zeros(N,N);
u = zeros(N,N);
v = zeros(N,N);
%%% SOLVING LOOP SIMILAR TO GAUSS-SIEDEL METHOD
for iter = 1:maxIt
%BOUNDARY CONDITIONS
w(1:N,N) = -2*psi(1:N,N-1)/(h^2) - U_wall*2/h; % Top
w(1:N,1) = -2*psi(1:N,2) /(h^2); % Bottom
w(1,1:N) = -2*psi(2,1:N) /(h^2); % Left
w(N,1:N) = -2*psi(N-1,1:N)/(h^2); % Right
%%% SINGLE STEP OF SOLVING VORTICITY TRANSPORT EQUATION
wp = w;
w(i,j) = wp(i,j) + (-1*(psi(i,j+1)-psi(i,j-1))/(2*h) .* (wp(i+1,i)-wp(i-1,j))/(2*h)+(psi(i+1,j)-psi(i-1,j))/(2*h) .* (wp(i,j+1)-wp(i,j-1))/(2*h)+...
mu/rho*(wp(i+1,j)+wp(i-1,j)-4*wp(i,j)+wp(i,j+1)+wp(i,j-1))/(h^2))*dt ;
%%% SINGLE STEP OF SOLVING VORTICITY - STREAM FUNCTION EQUATION
psi(i,j) = (w(i,j)*h^2 + psi(i+1,j) + psi(i,j+1) + psi(i,j-1) + psi(i-1,j))/4;
% CHECKING FOR CONVERGENCE

if iter > 40
    error = max(max(w - wp)) ;
    if error < maxe
        break;
    end
end
end
%%% CREATE VELOCITY FROM STREAM FUNCTION
u(2:N-1,N) = U_wall;
u(i,j) = (psi(i,j+1)-psi(i,j-1))/(2*h); % X-COMPONENT OF VELOCITY
v(i,j) = (-psi(i+1,j)+psi(i-1,j))/(2*h); %Y-COMPONENT OF VELOCITY

%%% PLOTS
x = 0:h:l; 
y = 0:h:l;
cm = hsv(ceil(100/0.7)); cm = flipud(cm(1:100,:));
%Plot of u
figure(1);
contourf(x,y,u',23,'LineColor','none');
title('U-velocity'); xlabel('x-location'); ylabel('y-location')
axis('equal',[0 l 0 l]); colormap(cm); colorbar('westoutside');
N = 1000;  xstart = max(x)*rand(N,1);  ystart = max(y)*rand(N,1);
[X,Y] = meshgrid(x,y);
figure(3);  h=streamline(X,Y,u',v',xstart,ystart,[0.1, 200]);
title('Stream Function');  xlabel('x-location');  ylabel('y-location')
axis('equal',[0 l 0 l]);  set(h,'color','k')
figure(4);
%plot of v 
contourf(x,y,v',23,'LineColor','none');
title('V-velocity'); xlabel('x-location'); ylabel('y-location')
axis('equal',[0 l 0 l]); colormap(cm); colorbar('westoutside');