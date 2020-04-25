clear all;clc;close all;
global G m
record=1; % shall we record video?


if record==1
    vid = VideoWriter('orbits_lepfrog_3D.mp4','MPEG-4');
    open(vid);
end

%% constants
N=16; % number of objects
G=1; % gravitational constant
m=ones(N,1)*1; % masses of the objects
max_iter_to_plot=10000; % how many points to trim the trajectory plot to

%% Initial conditions
% arranging objects -- optional, you can arrange them as you like
R=3; % radius of the circumference object will sit on initially

i=1:N; 
thetas=(360*i/N)';   % polar angles fo each object's initial position

X = R*cosd(thetas);         % Cartesian coordinates of the objects
Y = R*sind(thetas) ;        
Z = 0*thetas;               
        
r=[X Y Z];                  % position vectors of all masses

% Finding (configuration) Potential energy of each mass in the gravitational field due to all other masses
PE=zeros(N,1);
for i=1:N
    for j=1:N
        if i~=j
            PE(i)=PE(i)+G*m(i)./norm(r(j,:)-r(i,:)); % P.E. of mass i due to all j's
        end
    end
end
v0=sqrt(PE/2); % orbital Keplerian speeds for circular orbit (PE=2*KE) for each mass

% velocity components
vx = -v0.*sind(thetas);
vy = v0.*cosd(thetas);
vz = v0*0; % start motion in one plane
v=[vx vy vz];

%%%%%%% end initial conditions %%%%%%%%%%%


%% COM (barycentre)

COM=sum(r(1:N,1:3).*m(1:N),1)/sum(m); % location of barycentre

%% setup plotting
myplot=figure('Position',[100 100 850 850]);hold on; % figure window
COM_plot=plot3(COM(:,1),COM(:,2),COM(:,3),'*'); hold on;
axis([-10 10 -10 10 -10 10]/1.5);
% axis auto
xlabel('X');ylabel('Y');zlabel('Z');
daspect([1 1 1]); 
view(2); % change to view(3) for 3D projection
n_dts=360*128; % how many of timesteps per one period of the orbit

dt=max(2*pi*R/v0)/n_dts % timestep: n_dts steps per period of period

% plots for each orbit of each mass
for i=1:N
   body_traj(i)=plot3(X(i), Y(i),Z(i),'-'); hold on       % orbit line
   body(i)=plot3(X(i), Y(i), Z(i),'o', 'MarkerSize',8);    % mass marker
   set(body(i), 'MarkerFaceColor', body(i).Color, 'MarkerEdgeColor', body(i).Color)
end

t=0; % abs time
TITLE=title(['N=' num2str(N) ';  Gravitational N-body Simulation, Leapfrog DKD scheme']);
whitebg('black');

text(0,0,'COM');
E_txt=text(-5,-6,['Total Energy = ' ]);
time=text(2.5,3.5,['t = ' num2str(t)]);
grid on;
drawnow   


n=0; % iteration count

dt2=dt/2;% half timestep for leapfrog

while t<=60%isgraphics(myplot) % infinite time loop (while plot is open)
i=1:N;
r(i,1:3)=r(i,1:3)+v(i,1:3)*dt2; % half-step approximation

%% forces    
for i=1:N
   F(i,1:3)=zeros(1,3);
   for j=1:N
      if j~=i
         dr(i,:)=r(i,:)-r(j,:);
         r_cube=norm(dr(i,:)).^3; % scalar separation distance cubed
         F(i,1)=F(i,1)+G*m(i)*m(j)*dr(i,1)./r_cube; % x-comp of Force
         F(i,2)=F(i,2)+G*m(i)*m(j)*dr(i,2)./r_cube; % y-comp of Force
      end %if ends
   end % j
end % i; 
%% All net forces on all particles are calculated, find a v r Leapfrog DKD
i=1:N;
a(i,1:3)=-F(i,1:3)./m(i);
v(i,1:3)=v(i,1:3)+a(i,1:3)*dt;
r(i,1:3)=r(i,1:3)+v(i,1:3)*dt2;  
COM(1,1:3)=sum(r(1:N,:).*m(1:N),1)/sum(m,1); % new pos of COM

n=n+1; % number of total ode runs performed
t=t+dt;

% redrawind all plots
if mod(n,100)==0
    for ii=1:N
            body_traj(ii).XData=[body_traj(ii).XData r(ii,1)];
            body_traj(ii).YData=[body_traj(ii).YData r(ii,2)]; %    
            body_traj(ii).ZData=[body_traj(ii).ZData r(ii,3)]; % update plots    
            set(body(ii), 'XData',r(ii,1),'YData',r(ii,2),'ZData',r(ii,3)); % update mass' markers
            
        if length(body_traj(ii).XData)>=max_iter_to_plot % trimming the trajectory line to max_iter_to_plot
            body_traj(ii).XData=body_traj(ii).XData(end-max_iter_to_plot+1:end);
            body_traj(ii).YData=body_traj(ii).YData(end-max_iter_to_plot+1:end);
            body_traj(ii).ZData=body_traj(ii).ZData(end-max_iter_to_plot+1:end);
        end
    end  %%%%%%%%%%%%%%%%% loop ii   ends

COM(1:3)=sum(r(i,1:3).*m(i),1)/sum(m); % new location of barycentre
set(COM_plot,'XData', COM(1),'YData', COM(2),'ZData', COM(3)); 

t=t+dt; 
time.String=['t = ' num2str(round(t,3))];
drawnow

if record==1
    A=getframe(myplot);
    writeVideo(vid,A);
end


%% energy
% Kinetic Energy
V0=sqrt(sum(v(i,1:3).^2,2));
KE=sum(1/2*m.*V0.^2,1); % total Kintetic Energy


% Potential Energy
PE=zeros(N,1);
for i=1:N
    for j=1:N
        if i~=j
            PE(i)=PE(i)-G*m(i)./norm(r(j,:)-r(i,:)); % P.E. of mass i due to all j's
        end
    end
end
PE_tot=sum(PE);
E_txt.String=['Total Energy = ' num2str(KE+PE_tot)];
end
 end % while time loop ends

if record==1 
    close(vid); % close video file
end 
