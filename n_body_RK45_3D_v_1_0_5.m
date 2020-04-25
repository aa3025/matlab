clear all;clc;close all;
global G m
record=0; % shall we record video?


if record==1
    vid = VideoWriter('orbits_rk4_3D.mp4','MPEG-4');
    open(vid);
end

%% constants
N=16; % number of objects
G=1; % gravitational constant
m=ones(N,1)*1; % masses of the objects
n_dts_per_ode = 100; % number of timesteps to solve with ode solver during one iteration of while time loop
max_iter_to_plot=100000; % how many points to trim the trajectory plot to
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

inits=[X Y Z vx vy vz]; % initial conditions for each mass -- all in 1 column, 6 rows for ode solvers

%%%%%%% end initial conditions %%%%%%%%%%%


%% COM (barycentre)

COM=sum(r(1:N,1:3).*m(1:N),1)/sum(m); % location of barycentre

%% setup plotting
myplot=figure('Position',[100 100 850 850]);hold on; % figure window
COM_plot=plot3(COM(:,1),COM(:,2),COM(:,3),'*'); hold on;
axis([-10 10 -10 10 -10 10]*1);
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
end

t=0; % abs time
TITLE=title(['N=' num2str(N) ';  Gravitational N-body Simulation']);
whitebg('black');

text(0,0,'COM');
E_txt=text(0,-8,['Total Energy = ' ]);
time=text(2.5,3.5,['t = ' num2str(t)]);
grid on;
drawnow   

maxtime=dt*n_dts_per_ode; % time span for one ode45 call

time_span=0:dt:maxtime;
options = odeset('reltol',1e-13,'abstol',1e-12); % tolerance for ode solver


n=0; % iteration count


while isgraphics(myplot) % infinite DKD leapfrog time loop (while plot open)
    x_history=[];
    [TTT,XXX] = ode45(@nbody,time_span,inits,options);
    steps=length(TTT);
    x_history=XXX(:,:);

% get new coordinates from solution    
    i=1:N;
    XX(:,i)=x_history(:,i);
    YY(:,i)=x_history(:,i+N);
    ZZ(:,i)=x_history(:,i+2*N);
    
    VX(:,i)= x_history(:,i+3*N);
    VY(:,i)= x_history(:,i+4*N);
    VZ(:,i)= x_history(:,i+5*N);
    
    n=n+1; % number of total ode runs performed
    t=t+dt;

% redrawind all plots
    for ii=1:N
            body_traj(ii).XData=[body_traj(ii).XData XX(:,ii)'];
            body_traj(ii).YData=[body_traj(ii).YData YY(:,ii)']; %    
            body_traj(ii).ZData=[body_traj(ii).ZData ZZ(:,ii)']; % update plots    
            set(body(ii), 'XData',XX(end,ii),'YData',YY(end,ii),'ZData',ZZ(end,ii)); % update mass' markers
            
        if n*steps>=max_iter_to_plot % trimming the trajectory line to max_iter_to_plot
            body_traj(ii).XData=body_traj(ii).XData(end-max_iter_to_plot+1:end);
            body_traj(ii).YData=body_traj(ii).YData(end-max_iter_to_plot+1:end);
            body_traj(ii).ZData=body_traj(ii).ZData(end-max_iter_to_plot+1:end);
        end
    end  %%%%%%%%%%%%%%%%% loop ii   ends
  
r=[XX(end,:)' YY(end,:)' ZZ(end,:)'];
COM=sum(r(i,1:3).*m(i),1)/sum(m); % new location of barycentre
set(COM_plot,'XData', COM(end,1),'YData', COM(end,2),'ZData', COM(end,3)); 

inits=XXX(end,:)'; % use current state as initial conditions for next portion of orbit

t=t+maxtime-dt; 
time.String=['t = ' num2str(round(t,3))];
drawnow

if record==1
    A=getframe(myplot);
    writeVideo(vid,A);
end


%% energy
% Kinetic Energy
V0=sqrt(VX(end,:).^2+VY(end,:).^2+VZ(end,:).^2);
KE=sum(1/2*m'.*V0.^2,2); % total Kintetic Energy


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

 end % while time loop ends

if record==1 
    close(vid); % close video file
end 
