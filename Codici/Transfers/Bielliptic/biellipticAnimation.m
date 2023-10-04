% This script generates the animation for the bielliptic transfer between
% initial and final orbit.
% The maneuver uses two orbits to reach the ending point:
% - The first orbit enlarges the initial one so that the apocenter radius
%   is equal to r_a_m
% - The second orbit connects the point at r_a_m radius with the ending
%   point.
% This script also calculates the fuel mass ratio needed for the whole transfer 
% using the Tsiolkovsky equation.
% Parameters:
%   - Export (line 58):
%           -> 0: Plot animation without saving frames
%           -> 1: Save all frame as .png images
%   - Maximize (line 57);
%           -> 0: Open figure window with default size
%           -> 1: Open figure window maximized
%
% NOTE:
% The values used as apocenter radius are calculated in biellipticTests.m
% The file contains the different tests made and the results gathered,
% along with a brief analysis of the results obtained.
% As can be seen from that file, the optimal transfer is achieved when 
% the bielliptic maneuver starts at initial orbit's perigee and arrives at the ending point.

%% Animation engine

clear, clc
close all

addpath("..\..\Plot\","..\..\Utils\");

% -------------- ORBITS ---------------

% initial orbit
rr_i = [-1788.3462 -9922.9190 -1645.8335];
vv_i = [5.6510 -1.1520 -1.8710];
mu = 398600;
[a_i, e_i, i_i, OM_i, om_i, theta_i] = rv2parorb(rr_i, vv_i, mu);

% final orbit
a_f = 13290;
e_f = 0.3855;
i_f = Angle(0.9526);
OM_f = Angle(2.5510);
om_f = Angle(2.2540);
theta_f = Angle(3.0360);
[R_f, V_f] = parorb2rv(a_f, e_f, i_f, OM_f, om_f, theta_f, mu);

% -------------------------------------


% ------------- VARIABLES -------------

T_size = 100;                       % Amount of points per orbit
I_sp = 348;                         % Engine specific impulse (Engine: Merlin 1D+ Vacuum)

maximize = 1;
export = 0;
filename = '.\Export\Bielliptic';   % Filepath of saved frames (folder + frame name)
frame_duration = 0.0333;            % Number of seconds for each frame (in real time)
number_of_frames = 20;              % Number of frames for maneuver's animation

% -------------------------------------



% --------------- SETUP ---------------

if export
    export = answerFunction(export);
end
[theta_vect] = calculateThetaVect(mu, a_i, e_i, T_size);

if maximize
    figure("WindowState", "maximized");
else
    figure();
end
view(140,20);
hold on, grid on, axis equal
xlabel('X [Km]',FontSize=15);
ylabel('Y [Km]',FontSize=15);
zlabel('Z [Km]',FontSize=15);

% Earth plot
s = earthPlot;

one_orbit_time = frame_duration*T_size;
orbit_duration = 2*pi*sqrt((a_i^3)/mu);

factor = one_orbit_time/orbit_duration;

earth_rotation_duration = 24*60*60;
earth_rotation_duration = earth_rotation_duration*factor;
angular_vel = 360/earth_rotation_duration;
frame_rotation_angle = angular_vel*frame_duration;

% Filename
name = append(filename,"_1.png");

% Result variables
dv_tot = 0;
dt_tot = 0;

% Global frame counter
j = 0; 

% Apocenter radius: where the maneuver will be executed
r_a_m = 57000;

% -------------------------------------



% -------------- ENGINE ---------------

% initial and final orbits
theta_vect1=linspace(0,2*pi,500);
[rr_i_vect, ~] = parorb2rv(a_i, e_i, i_i, OM_i, om_i, theta_vect1, mu);
[rr_f_vect, ~] = parorb2rv(a_f, e_f, i_f, OM_f, om_f, theta_vect1, mu);
plot3(rr_i_vect(:,1),rr_i_vect(:,2),rr_i_vect(:,3),Color='#ff9500',LineWidth=0.4);
plot3(rr_f_vect(:,1),rr_f_vect(:,2),rr_f_vect(:,3),Color='#ff9500',LineWidth=0.4);

plot3(rr_i(1),rr_i(2),rr_i(3),'xr',LineWidth=2);
[rr_f_point, ~] = parorb2rv(a_f, e_f, i_f, OM_f, om_f, theta_f, mu);
plot3(rr_f_point(1),rr_f_point(2),rr_f_point(3),'xr',LineWidth=2);



% Framing
max_frame=r_a_m;

set(gca,'XLim',[-max_frame max_frame])
set(gca,'YLim',[-max_frame max_frame])
set(gca,'ZLim',[-max_frame max_frame])


% Calculate first elliptic transfer orbit
[R_p_i, V_p_i] = parorb2rv(a_i, e_i, i_i, OM_i, om_i, Angle(0), mu);
a_tras_1 = (norm(R_p_i) + r_a_m)/2;
e_tras_1 = (r_a_m - norm(R_p_i))/(r_a_m + norm(R_p_i));


% Plot along initial orbit to reach first maneuver point
angle1 = theta_i;
angle2 = Angle(0);
% First maneuver: Initial orbit widening
[rr_vect_complete, ~] = parorb2rv(a_i,e_i,i_i,OM_i,om_i,theta_vect,mu); % complete orbit
[rr_semiVect] = semiOrb(rr_vect_complete,angle1,angle2,a_i, e_i, i_i, OM_i, om_i, mu);

type=0;
plot_temp=plot3(rr_vect_complete(:,1),rr_vect_complete(:,2),rr_vect_complete(:,3),'--m');
for k=1:size(rr_semiVect,1)
    j=j+1;
    rr_orb=rr_semiVect(1:k,:);
    exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle,"blu");
%     exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle,"#bf00ff");
end
delete(plot_temp)



% Animate orbit change

[R_p_tras_1, V_p_tras_1] = parorb2rv(a_tras_1, e_tras_1, i_i, OM_i, om_i, angle2, mu);
dv_1 = norm(V_p_tras_1 - V_p_i);

dv_1_step=dv_1./number_of_frames;
dv_1_vect=(1:number_of_frames)'.*dv_1_step;
dv_1_vect = flip(dv_1_vect);

type=1;
plot_temp=plot3(rr_semiVect(end,1),rr_semiVect(end,2),rr_semiVect(end,3),'or');
for k=2:number_of_frames
    j=j+1;
    vv_1_i_post=V_p_tras_1+dv_1_vect(k, :);
    [a_tras, e_tras, i_tras, OM_tras, om_tras, theta_tras] = rv2parorb(R_p_tras_1, vv_1_i_post, mu);
    [T_size_new] = sizeOfTime(T_size,orbit_duration,a_tras);
    [theta_vect] = calculateThetaVect(mu, a_tras, e_tras, T_size_new);
    [rr_peri_tras, ~] = parorb2rv(a_tras, e_tras, i_tras, OM_tras, om_tras, theta_vect, mu);
    exportPng(rr_peri_tras,j,filename,export,type,s,frame_rotation_angle);
end
delete(plot_temp)
scatter3(rr_semiVect(end,1),rr_semiVect(end,2),rr_semiVect(end,3), 'blue', 'filled'); % Plot maneuver point


% Update results
dt_1 = timeCalc(a_i, e_i, mu, angle1, angle2);
dv_tot = dv_tot + dv_1;
dt_tot = dt_tot + dt_1;



% Second maneuver: Transfer orbit to reach ending point
angle1=angle2;          % relative to the 1st orbit, rad
angle2=angle1+pi;       % relative to the 1st orbit, rad

% Plot orbit to reach second maneuver point
[T_size_new] = sizeOfTime(T_size,orbit_duration,a_tras_1);
[theta_vect] = calculateThetaVect(mu, a_tras_1, e_tras_1, T_size_new);
[rr_vect_complete, ~] = parorb2rv(a_tras_1,e_tras_1,i_i,OM_i,om_i,theta_vect,mu); % complete orbit
[rr_vect_complete] = VectorCorrector(rr_vect_complete,3000);
[rr_semiVect] = semiOrb(rr_vect_complete,angle1,angle2,a_tras, e_tras, i_i, OM_i, om_i, mu);

type=0;
plot_temp=plot3(rr_vect_complete(:,1),rr_vect_complete(:,2),rr_vect_complete(:,3),'--m');
for k=1:size(rr_semiVect,1)
    j=j+1;
    rr_orb=rr_semiVect(1:k,:);
    exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle,"blu");
%     exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle,"#0ba600");
end
delete(plot_temp)



% Calculate second transfer orbit
[R_a_m_1, V_a_m_1] = parorb2rv(a_tras_1, e_tras_1, i_i, OM_i, om_i, Angle(pi), mu);
[a_tras_2, e_tras_2, i_tras_2, OM_tras_2, om_tras_2, th_1, th_2] = twoPointsOrbit(R_a_m_1, R_f, 0);
[R_a_m_2, V_a_m_2] = parorb2rv(a_tras_2, e_tras_2, i_tras_2, OM_tras_2, om_tras_2, th_1, mu);


% Animate orbit change
[i_vect,OM_vect,om_vect]=parametersSpread(number_of_frames,i_i,i_tras_2,OM_i,OM_tras_2,om_i,om_tras_2);
[a_vect,e_vect,~]=parametersSpread(number_of_frames,a_tras_1,a_tras_2,e_tras_1,e_tras_2,0,0);

type=1;
plot_temp=plot3(R_a_m_1(1),R_a_m_1(2),R_a_m_1(3),'or');
for k=2:number_of_frames % without initial orbit
    j=j+1;  % global
    [T_size_new] = sizeOfTime(T_size,orbit_duration,a_vect(k));
    [theta_vect] = calculateThetaVect(mu, a_vect(k), e_vect(k), T_size_new);
    [rr_orb, ~] = parorb2rv(a_vect(k),e_vect(k),i_vect.value(k),OM_vect.value(k),om_vect.value(k),theta_vect,mu);
    [rr_orb] = VectorCorrector(rr_orb,4000);
    exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle);
end
delete(plot_temp)
scatter3(rr_semiVect(end,1),rr_semiVect(end,2),rr_semiVect(end,3), 'blue', 'filled'); % Plot maneuver point


% Update results
dv_2 = norm(V_a_m_2 - V_a_m_1);
dt_2 = timeCalc(a_tras_1, e_tras_1, mu, angle1, angle2);
dv_tot = dv_tot + dv_2;
dt_tot=dt_tot+dt_2;


% Plot along orbit to reach ending point
angle1=angle2;          
angle2=th_2;   

[T_size_new] = sizeOfTime(T_size,orbit_duration,a_tras_2);
[theta_vect] = calculateThetaVect(mu, a_tras_2, e_tras_2, T_size_new);
[rr_vect_complete, ~] = parorb2rv(a_tras_2,e_tras_2,i_tras_2,OM_tras_2,om_tras_2,theta_vect,mu); % complete orbit
[rr_semiVect] = semiOrb(rr_vect_complete,angle1,angle2,a_tras_2,e_tras_2,i_tras_2,OM_tras_2,om_tras_2, mu);

type=0;
plot_temp=plot3(rr_vect_complete(:,1),rr_vect_complete(:,2),rr_vect_complete(:,3),'--m');
for k=1:size(rr_semiVect,1)
    j=j+1;
    rr_orb=rr_semiVect(1:k,:);
    exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle,"blu");
%     exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle,"#ff5500");
end
delete(plot_temp)


% Animate orbit change
[R_m_2, V_m_2] = parorb2rv(a_tras_2, e_tras_2, i_tras_2, OM_tras_2, om_tras_2, th_2, mu);
[i_vect,OM_vect,om_vect]=parametersSpread(number_of_frames,i_tras_2,i_f,OM_tras_2, ...
    OM_f,om_tras_2,om_f);
[a_vect,e_vect,~]=parametersSpread(number_of_frames,a_tras_2,a_f,e_tras_2,e_f,0,0);

type=1;
plot_temp=plot3(R_m_2(1),R_m_2(2),R_m_2(3),'or');
for k=2:number_of_frames % without initial orbit
    j=j+1;  % global
    [T_size_new] = sizeOfTime(T_size,orbit_duration,a_vect(k));
    [theta_vect] = calculateThetaVect(mu, a_vect(k), e_vect(k), T_size_new);
    [rr_orb, ~] = parorb2rv(a_vect(k),e_vect(k),i_vect.value(k),OM_vect.value(k),om_vect.value(k),theta_vect,mu);
    exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle);
end
delete(plot_temp)
scatter3(rr_semiVect(end,1),rr_semiVect(end,2),rr_semiVect(end,3), 'blue', 'filled'); % Plot maneuver point



% Update results
dv_3 = norm(V_f - V_m_2);
dt_3 = timeCalc(a_tras_2, e_tras_2, mu, angle1, angle2);
dv_tot = dv_tot + dv_3;
dt_tot=dt_tot+dt_3;



if export
    warning("Done!")
end

% Display results
disp("Delta V tot = "+dv_tot+" km/s")
fprintf("\t--> delta V1 = "+norm(dv_1)+" km/s\n")
fprintf("\t--> delta V2 = "+norm(dv_2)+" km/s\n")
fprintf("\t--> delta V3 = "+norm(dv_3)+" km/s\n")

disp("Delta t tot = "+dt_tot+" s");
fprintf("\t--> delta t1 = "+dt_1+" s\n")
fprintf("\t--> delta t2 = "+dt_2+" s\n")
fprintf("\t--> delta t3 = "+dt_3+" s\n")

% -------------------------------------



% ------------ TSIOLKOVSKY ------------

g = 9.807;

% dv must be converted to m/s!
dv_tot_ms = dv_tot*1000;

% Calculate mass ratio from ln(m_ratio)=dv/(I_sp*g)
m_ratio=exp(dv_tot_ms/(I_sp*g));

% Calculate fuel mass

% m_ratio=m0/mf;
% m0=mf+mc;             where mf is dry mass, mc is fuel mass and m0 is wet mass
% m_ratio=(mf+mc)/mf;
% m_ratio=1+mc/mf;
% m_ratio-1=mc/mf;

mc = m_ratio-1;     % mc=mc(mf)

fprintf("\nPropellant mass = ("+mc+" kg)*dryMass, with a single stage ("+I_sp+" s impulse) engine")

% fuelMass=mc*dryMass
% totalMass=mc*dryMass+dryMass=dryMass(mc+1)

% %: mc*dryMass/dryMass(mc+1)=mc/(mc+1)

fprintf("\nMass percentage (fuelMass/totalMass)*100 = "+mc/(mc+1)*100+"%%\n")

% -------------------------------------



% ---------- LOCAL FUNCTIONS ----------

function exportPng(rr_vect,k,filename,export,type,s,frame_rotation_angle,color)
    % Generate plot and save it as png, if export is True
    
    x = rr_vect(:,1);
    y = rr_vect(:,2);
    z = rr_vect(:,3);
    
    if ~type
        plot3(x,y,z,'Color',color,LineWidth=1.3);
        frame1 = plot3(x(end),y(end),z(end),'or');
        rotate(s, [0 0 1],frame_rotation_angle,[0 0 0])
    else
        frame1 = plot3(x,y,z,'--m');
    end

    pause(0.1)
    % Save frame as .png file
    if export
        name = append(filename,"_",num2str(k),".png");
        exportgraphics(gcf,name,'Resolution',300)
    end
    
    delete(frame1)
end


function [param_1,param_2,param_3] = parametersSpread(num,param_1_start,param_1_end,param_2_start,param_2_end,param_3_start,param_3_end)
    % Return vectors from param_start to param_end with num elements
    param_1 = linspace(param_1_start,param_1_end,num);
    param_2 = linspace(param_2_start,param_2_end,num);
    param_3 = linspace(param_3_start,param_3_end,num);
end


function [export] = answerFunction(export)
    answer = 0;
    while ~answer
        reply = input('Are you sure to create PNGs? [Y/N]:','s');
        if isempty(reply)
            reply = 'Y';
        end
        if strlength(reply) == 1
            if reply == 'n' || reply == 'N'
                warning("'export' set to 0")
                export = 0;
                answer = 1;
            elseif reply == 'y' || reply == 'Y'
                warning("Generating PNGs...")
                answer = 1;
            end
        end
    end
end

function [T_size_new] = sizeOfTime(T_size,orbit_duration,a)
% relative T_size for realistic orbital times
    mu=398600;
    orbit = 2*pi*sqrt((a^3)/mu);
    ratio=orbit/orbit_duration;
    T_size_new=ceil(T_size*ratio);
end

function [rr] = VectorCorrector(rr,maxDistance)
% fix for maneuver animation (graphic issue)
numberOfPoints=size(rr,1);
    for i=1:numberOfPoints-2
        distance=norm(rr(i+1,:)-rr(i,:));
        if distance>maxDistance
            rr(i+1,:)=rr(i,:);
        end
    end
end

% -------------------------------------