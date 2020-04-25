function [dx] = nbody(t,x)
global G m

N = length(m); % The number of bodies
dx = zeros(6*N,1); % for each particle we need 6 equations: 3 for position, 3 for velocity; 

% computing distances r between bodies
X=x(1:N);
Y=x(N+1:2*N);
Z=x(2*N+1:3*N);

DX=distance(X); % external function below
DY=distance(Y);
DZ=distance(Z);

r=(DX.^2+DY.^2+DZ.^2).^0.5; % matrix of distances between all particles

%% accelerations
   M=m*ones(1,N); % square matrix of masses, each row is m(i) repeated along columns (for next line)
   a=-G*M./r.^2; % accelerations of each mass
   
   % get rid of diagonal entries (Inf - forces on bodies on itself)
   a(1:N+1:N*N) = 0;
    
   ax =a.* (DX./ r);
   ay =a.* (DY./ r);
   az =a.* (DZ./ r);
   
   ax(1:N+1:N*N) = 0;
   ay(1:N+1:N*N) = 0;
   az(1:N+1:N*N) = 0;
   
   % net accel on mass m(i) -- adding all a's along column number (direction 2)
   a_sx=sum(ax,2)';
   a_sy=sum(ay,2)';
   a_sz=sum(az,2)';

%% velocities

% Filling in the velocities from IC's
for i = 1:3*N;
    dx(i) = x(3*N+i);    % dx is [vx1 vx2...vxN vy1 vy2...vyN vz1 vz2...vzN ax1 ax2 ...axN ay...etc az...etc]
end

for i=1:N
    dx(3*N+i)=a_sx(i)';
    dx(4*N+i)=a_sy(i)';
    dx(5*N+i)=a_sz(i)';
end
end



function [dist]=distance(X)
% this function calculates distance between two coordinates in a given vector X
L=length(X);
dist=zeros([1 L]);

j=1:L;
for k=1:L;
  dist(j(k),j)=X(j(k))-X(j);
end
 
end