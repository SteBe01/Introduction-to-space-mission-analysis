% This script is used to generate the animation (fixed camera) with
% the satellite moving across the initial orbit

clear, clc
close all

export=0;
k=0;
filename = '.\Export\Camera';
addpath(".\Utils\",".\Plot\");

% initial orbit
rr_i = [-1788.3462 -9922.9190 -1645.8335];
vv_i = [5.6510 -1.1520 -1.8710];
mu = 398600;
[a_i, e_i, i_i, OM_i, om_i, theta_i] = rv2parorb(rr_i, vv_i, mu);

% final orbit
a_f = 13290;
e_f = 0.3855;
i_f = 0.9528;
OM_f = 2.5510;
om_f = 2.2540;
theta_f = 3.0360;

figure("WindowState", "maximized");
s=earthPlot;
hold on, axis off, axis equal
frame_duration=0.0333;
T_size=60;

one_orbit_time = frame_duration*T_size;
orbit_duration = 2*pi*sqrt((a_i^3)/mu);
factor = one_orbit_time/orbit_duration;
earth_rotation_duration = 24*60*60+56*60+04; % 1 day == 24h 56m 04s
earth_rotation_duration = earth_rotation_duration*factor;
angular_vel = 360/earth_rotation_duration;
frame_rotation_angle = angular_vel*frame_duration;

[theta_vect] = calculateThetaVect(mu, a_i, e_i, T_size);
[rr, vv] = parorb2rv(a_i, e_i, i_i, OM_i, om_i, theta_vect, mu);

h=cross(rr_i,vv_i);
h=h/norm(h);

% framing settings
[apo, ~] = parorb2rv(a_i, e_i, i_i, OM_i, om_i, pi, mu);
max_frame=2*norm(apo);
set(gca,'XLim',[-max_frame max_frame])
set(gca,'YLim',[-max_frame max_frame])
set(gca,'ZLim',[-max_frame max_frame])
camzoom(3)

view([rr(1,1),-rr(1,2),rr(1,3)])

for i=1:size(rr,1)-1
    k=k+1;
    plot=plot3(rr([1:i],1),rr([1:i],2),rr([1:i],3),Color='blue');
    plot1=plot3(rr(i,1),rr(i,2),rr(i,3),'or');

    angle=360*((theta_vect(i+1)-theta_vect(i))/(2*pi));
    camtarget([rr(i,1),rr(i,2),rr(i,3)])

    camorbit(angle,0,'data',h)

    pause(0.1)
    if export
        name = append(filename,"_",num2str(k),".png");
        % export framing
        a = annotation('rectangle',[0.4 0 0.6 0.8],'Color','#FFFFFE');
        exportgraphics(gcf,name,'Resolution',300)
        delete(a)
    end

    delete(plot)
    delete(plot1)
    rotate(s, [0 0 1],frame_rotation_angle,[0 0 0])
end