% We going to build a "snake", which will be just a set of dots(joints)
% connected with short lines. First dot ("head") will react to the mouse position, and all
% other dots will follow along the chain. Alex Pedcenko 2017, alex.pedcenko(at)coventry.ac.uk.
clc,close all
N=50; % number of "dots/joints" in snake
dl=20; % distance between joints
x=0:dl:dl*(N-1); % assigning initial positions for joints for x
y=zeros(1,N); %... and for y
r=[x' y']; % position vector for each joint (i.e. table of coords  , i.e. r = [x  y] )
figure('units','normalized','outerposition',[0 0 1 1]) % bring up the empty fullscreen plot window
joints=plot(r(:,1),r(:,2),'o','MarkerFaceColor','m', 'MarkerEdgeColor','m'); % plot joints of the snake with magenta
hold on
links=plot(r(:,1),r(:,2),'-m'); % plot lines between joints
daspect([1 1 1]);
% find min/max of coordinates of the snake
maxx=max(abs(r(:,1)))/4;
maxy=max(abs(r(:,2)))/4;
axis([-maxx*1.8 maxx*1.8 -maxx maxx]) % set min max on axes
axis manual % do not rescale axes
h=title('Hover mouse over the plot');
% set what happens on MouseMove event (disable below line to make mouse
% response to click only)
set(gcf,'WindowButtonMotionFcn', 'drawnow'); % 'drawnow' is executed on mouse move
% Main loop here for updating snake coords
while ~isempty(get(groot,'CurrentFigure')); % inifinite loop while figure window exists
C=get(gca,'CurrentPoint');% get mouse coordinates
r(1,:)=C(1,1:2); % Assign ccordinates of mouse to the "head" of the snake
for i=2:N % loop along joints from the 2nd to last
    dr=r(i,:)-r(i-1,:);
    r_roof=dr./norm(dr); % unit vector of the link between this particle and the previous one
    
    % update coordinates of the current joint as being "dl" units away from
    % the previous joint along the direction r_roof 
    
    r(i,1)=r(i-1,1)+r_roof(1)*dl; % x
    r(i,2)=r(i-1,2)+r_roof(2)*dl; % y
drawnow
end
 % update plots with new coordinates of joints and links
try % with "try... end" , if the figure is closed the next lines will not produce an error
    set(joints,'XData',r(:,1),'YData',r(:,2));
    set(links,'XData',r(:,1),'YData',r(:,2));
end
end
