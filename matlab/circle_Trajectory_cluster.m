%%
clear 
clc
clf

start = zeros(4,1);
state = zeros(4,1);
R = 50;
interR = 16;
aplha = linspace(0, 2*pi, 200);
rmin = -pi/2;
rmax = pi/2;
hmin = -pi/4;
hmax = pi/4;
step = (hmax-hmin)/4;
figure(1)
plot(R*cos([0:0.1:2*pi]),R*sin([0:0.1:2*pi]),'r--'), hold on ;
Xstate = [];
j = 0;

for i = rmin:step:rmax
%     for j = hmin:step:hmax
        state = round( [ R*cos(i); R*sin(i); i+j; 0 ] , 8);     % 取4位精度
        if abs( state(3)-start(3)) > pi
            state(3) = gui(state(3));
            start(3) = gui(start(3));
        end
        Xstate(:,end+1) = state;
        [ Path, S, status] = trajectory_generate(50, start,state);
        if status == 0
            continue;
        end
        plot( Path(1,:), Path(2,:) ), hold on;
%     end
end

start(3) = -pi;
for i = rmax:step:2*pi+rmin
%     for j = hmin:step:hmax
        state = round( [ R*cos(i); R*sin(i); i+j; 0 ] , 8);     % 取4位精度
        if abs( state(3)-start(3)) > pi
            state(3) = gui(state(3));
            start(3) = gui(start(3));
        end
        Xstate(:,end+1) = state;
        [ Path, S, status] = trajectory_generate(50, start,state);
        if status == 0
            continue;
        end
        plot( Path(1,:), Path(2,:) ), hold on;
%     end
end


function x = gui( x )
    if x<-pi
        x = x + (fix(x/(2*pi))+1.0)*2*pi;
    elseif x>=pi
        x = x - (fix(x/(2*pi))+1.0)*2*pi;
    end        
end