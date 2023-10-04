% This script generates the animation for the standard transfer between
% initial and final orbit, with optimization for dv or dt.
% It also calculates the fuel mass ratio needed for the whole transfer 
% using the Tsiolkovsky equation.
% Parameters:
%   - Optimization (line 45): 
%           -> 0: transfer with lowest dt
%           -> 1: transfer with lowest dv
%   - Export (line 50):
%           -> 0: Plot animation without saving frames
%           -> 1: Save all frame as .png images
%   - Maximize (line 49);
%           -> 0: Open figure window with default size
%           -> 1: Open figure window maximized
%

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

% -------------------------------------


% ------------- VARIABLES -------------

optimization = 0;                   
T_size = 300;                       % Amount of points per orbit
I_sp = 348;                         % Engine specific impulse (Engine: Merlin 1D+ Vacuum)

maximize = 1;
export = 0;
filename = '.\Export\Standard';     % Filepath of saved frames (folder + frame name)
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
view(115,20);
hold on, grid on, axis equal
xlabel('X [Km]',FontSize=15);
ylabel('Y [Km]',FontSize=15);
zlabel('Z [Km]',FontSize=15);

% Earth plot
s = earthPlot;

one_orbit_time = frame_duration*T_size;
orbit_duration = 2*pi*sqrt((a_i^3)/mu);

factor = one_orbit_time/orbit_duration;

earth_rotation_duration = 24*60*60+56*60+04; % 1 day == 24h 56m 04s
earth_rotation_duration = earth_rotation_duration*factor;
angular_vel = 360/earth_rotation_duration;
frame_rotation_angle = angular_vel*frame_duration;

% Framing
[apo_i, ~] = parorb2rv(a_i, e_i, i_i, OM_i, om_i, pi, mu);
[apo_f, ~] = parorb2rv(a_f, e_f, i_f, OM_f, om_f, pi, mu);
max_frame = max(norm(apo_i),norm(apo_f));

set(gca,'XLim',[-max_frame max_frame])
set(gca,'YLim',[-max_frame max_frame])
set(gca,'ZLim',[-max_frame max_frame])

% Filename
name = append(filename,"_1.png");

% Result variables
dv_tot = 0;
dt_tot = 0;

% Global frame counter
j = 0; 

% -------------------------------------



% -------------- ENGINE ---------------

% initial and final orbits
theta_vect1=linspace(0,2*pi,500);
[rr_i_vect, ~] = parorb2rv(a_i, e_i, i_i, OM_i, om_i, theta_vect1, mu);
[rr_f_vect, ~] = parorb2rv(a_f, e_f, i_f, OM_f, om_f, theta_vect1, mu);
plot3(rr_i_vect(:,1),rr_i_vect(:,2),rr_i_vect(:,3),Color='#ff9500',LineWidth=0.4);
plot3(rr_f_vect(:,1),rr_f_vect(:,2),rr_f_vect(:,3),Color='#ff9500',LineWidth=0.4);

plot3(rr_i(1),rr_i(2),rr_i(3),'xr',LineWidth=2,MarkerSize=6);
[rr_f_point, ~] = parorb2rv(a_f, e_f, i_f, OM_f, om_f, theta_f, mu);
plot3(rr_f_point(1),rr_f_point(2),rr_f_point(3),'xr',LineWidth=2,MarkerSize=6);



% First maneuver: Plane change
[dv_1_plane, dv_2_plane, theta_cambioPiano_1, theta_cambioPiano_2, om_1] = planeChange(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu);

% Find plane change angle based on optimization
time_forced = 0;        % If dv is equal in both points, optimize for dt
if optimization
    if dv_2_plane<dv_1_plane
        dv_1 = dv_2_plane;
        theta_cambioPiano = theta_cambioPiano_2;
    elseif dv_2_plane>dv_1_plane
        dv_1 = dv_1_plane;
        theta_cambioPiano = theta_cambioPiano_1;
    else
        time_forced = 1;
    end
end
if ~optimization || time_forced
    dummy_1 = theta_cambioPiano_1-theta_i;
    dummy_2 = theta_cambioPiano_2-theta_i;
    if dummy_1<pi
        theta_cambioPiano = theta_cambioPiano_1;
        dv_1 = dv_1_plane;
    else
        theta_cambioPiano = theta_cambioPiano_2;
        dv_1 = dv_2_plane;
    end
end


% Plot along orbit to reach first maneuver angle
angle1 = theta_i;                 % [rad], relative to the 1st orbit
angle2 = theta_cambioPiano;       % [rad], relative to the 1st orbit
dt_1 = timeCalc(a_i, e_i, mu, angle1, angle2);

[rr_vect_complete, ~] = parorb2rv(a_i,e_i,i_i,OM_i,om_i,theta_vect,mu);
[rr_semiVect] = semiOrb(rr_vect_complete,angle1,angle2,a_i, e_i, i_i, OM_i, om_i, mu);

type = 0;
plot_temp = plot3(rr_vect_complete(:,1),rr_vect_complete(:,2),rr_vect_complete(:,3),'--m');
for k = 1:size(rr_semiVect,1)
    j = j+1;
    rr_orb = rr_semiVect(1:k,:);
    exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle,"blu");
%     exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle,"#bf00ff");
end
delete(plot_temp)


% Animate orbit change
[i_vect,OM_vect,om_vect] = parametersSpread(number_of_frames,i_i,i_f,OM_i,OM_f,om_i,om_1);

type = 1;
plot_temp = plot3(rr_semiVect(end,1),rr_semiVect(end,2),rr_semiVect(end,3),'or');
for k = 2:number_of_frames 
    j = j+1;  
    [rr_orb, ~] = parorb2rv(a_i,e_i,i_vect.value(k),OM_vect.value(k),om_vect.value(k),theta_vect,mu);
    exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle);
end
delete(plot_temp)
scatter3(rr_semiVect(end,1),rr_semiVect(end,2),rr_semiVect(end,3), 'blue', 'filled'); % Plot maneuver point


% Update results
dv_tot = dv_tot+norm(dv_1);
dt_tot = dt_tot+dt_1;



% Second maneuver: Pericenter anomaly change
[dv_2, theta_2_a, theta_2_b, om_2] = omChange(a_i, e_i, om_1, om_f, mu);

% Since dv is equal in both point, find angle that optimizes dt
angle1 = theta_cambioPiano;
dummy = (2*pi-theta_2_b)-angle1;
dummy2 = (2*pi-theta_2_a)-angle1;
if dummy>dummy2
    angle2 = 2*pi-theta_2_a;
else
    angle2 = 2*pi-theta_2_b;
end

dt_2 = timeCalc(a_i, e_i, mu, angle1, angle2);


% Plot along orbit to reach second maneuver angle
[rr_vect_complete, ~] = parorb2rv(a_i,e_i,i_f,OM_f,om_1,theta_vect,mu); % complete orbit
[rr_semiVect] = semiOrb(rr_vect_complete,angle1,angle2,a_i, e_i, i_f, OM_f, om_1, mu);

type = 0;
plot_temp = plot3(rr_vect_complete(:,1),rr_vect_complete(:,2),rr_vect_complete(:,3),'--m');
for k = 1:size(rr_semiVect,1)
    j = j+1;
    rr_orb = rr_semiVect(1:k,:);
    exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle,"blu");
%     exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle,"#0000ff");
end
delete(plot_temp)


% Animate orbit change
[i_vect,OM_vect,om_vect] = parametersSpread(number_of_frames,0,0,0,0,om_1,om_2);

type = 1;
plot_temp = plot3(rr_semiVect(end,1),rr_semiVect(end,2),rr_semiVect(end,3),'or');
for k = 2:number_of_frames 
    j = j+1;
    [rr_orb, ~] = parorb2rv(a_i,e_i,i_f,OM_f,om_vect.value(k),theta_vect,mu);
    exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle);
end
delete(plot_temp)
scatter3(rr_semiVect(end,1),rr_semiVect(end,2),rr_semiVect(end,3), 'blu', 'filled'); % Plot maneuver point

% Update results
dv_tot = dv_tot+norm(dv_2);
dt_tot = dt_tot+dt_2;



% Third maneuver: First bitangent burn to change shape
angle1 = 2*pi-angle2;

% Find best solution based on optimization
[dv_3_1, dv_4_1, dt_4_1] = orbitShapeChange(a_i, e_i, a_f, e_f, om_2, om_f, mu, OM_f, i_f, 1);
[dv_3_2, dv_4_2, dt_4_2] = orbitShapeChange(a_i, e_i, a_f, e_f, om_2, om_f, mu, OM_f, i_f, 0);
dv_change_ae_1 = norm(dv_3_1)+norm(dv_4_1);
dv_change_ae_2 = norm(dv_3_2)+norm(dv_4_2);

time_forced = 0;
if optimization
    if dv_change_ae_1<dv_change_ae_2
        choice = 1;
        angle2 = Angle(0);
    elseif dv_change_ae_1>dv_change_ae_2
        choice = 0;
        angle2 = Angle(pi);
    else
        time_forced = 1;
    end
end
if ~optimization || time_forced
    if dt_4_1<dt_4_2
        choice = 1;
        angle2 = Angle(0);
    else
        choice = 0;
        angle2 = Angle(pi);
    end
end

dt_3 = timeCalc(a_i, e_i, mu, angle1, angle2);


% Plot along orbit to reach third maneuver angle
[rr_vect_complete, ~] = parorb2rv(a_i,e_i,i_f,OM_f,om_2,theta_vect,mu); % complete orbit
[rr_semiVect] = semiOrb(rr_vect_complete,angle1,angle2,a_i, e_i, i_f, OM_f, om_2, mu);

type = 0;
plot_temp = plot3(rr_vect_complete(:,1),rr_vect_complete(:,2),rr_vect_complete(:,3),'--m');
for k = 1:size(rr_semiVect,1)
    j = j+1;
    rr_orb = rr_semiVect(1:k,:);
    exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle,"blu");
%     exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle,"#00b3a4");
end
delete(plot_temp)


% Animate orbit change
[dv_3, dv_4, dt_4] = orbitShapeChange(a_i, e_i, a_f, e_f, om_2, om_f, mu, OM_f, i_f, choice);

if choice
    [rr_1_i, vv_1_i] = parorb2rv(a_i, e_i, i_f, OM_f, om_2, 0, mu); % 1 = pericentro, 2 = apocentro
else
    [rr_1_i, vv_1_i] = parorb2rv(a_i, e_i, i_f, OM_f, om_2, pi, mu); % 1 = apocentro, 2 = pericentro
end

dv_3_step = dv_3./number_of_frames;
dv_3_vect = (1:number_of_frames)'.*dv_3_step;

type = 1;
plot_temp = plot3(rr_semiVect(end,1),rr_semiVect(end,2),rr_semiVect(end,3),'or'); 
for k = 2:number_of_frames 
    j = j+1;
    vv_1_i_post = vv_1_i+dv_3_vect(k, :);
    [a_tras, e_tras, i_tras, OM_tras, om_tras, theta_tras] = rv2parorb(rr_1_i, vv_1_i_post, mu);
    [T_size_new] = sizeOfTime(T_size,orbit_duration,a_tras);
    [theta_vect] = calculateThetaVect(mu, a_tras, e_tras, T_size_new);
    [rr_peri_tras, ~] = parorb2rv(a_tras, e_tras, i_tras, OM_tras, om_tras, theta_vect, mu);
    exportPng(rr_peri_tras,j,filename,export,type,s,frame_rotation_angle);
end
delete(plot_temp)
scatter3(rr_semiVect(end,1),rr_semiVect(end,2),rr_semiVect(end,3), 'blue', 'filled'); % Plot maneuver point


% Update results
dv_tot = dv_tot+norm(dv_3);
dt_tot = dt_tot+dt_3;
dt_tot = dt_tot+dt_4;



% Fourth maneuver: Second bitangent burn to change shape
[p1,p2,~,angle1] = intersection(a_i,e_i,i_f,OM_f,om_2,a_tras,e_tras,i_tras,OM_tras,om_tras,mu,1e-6,angle2,0);
if angle1-pi<angle1
    angle1 = Angle(pi);
else
    angle1 = Angle(0);
end
angle2 = angle1+pi;


% Plot along orbit to reach fourth maneuver angle
[rr_vect_complete, ~] = parorb2rv(a_tras,e_tras,i_f,OM_f,om_f,theta_vect,mu); % complete orbit
[rr_semiVect] = semiOrb(rr_vect_complete,angle1,angle2,a_tras, e_tras, i_f, OM_f, om_f, mu);

type = 0;
plot_temp = plot3(rr_vect_complete(:,1),rr_vect_complete(:,2),rr_vect_complete(:,3),'--m');
for k = 1:size(rr_semiVect,1)
    j = j+1;
    rr_orb = rr_semiVect(1:k,:);
    exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle,"blu");
%     exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle,"#0ba600");
end
delete(plot_temp)


% Animate orbit change
angle = theta_tras+pi;
[rr_2_tras, vv_2_tras] = parorb2rv(a_tras, e_tras, i_tras, OM_tras, om_tras, angle, mu);

dv_4_step = dv_4./number_of_frames;
dv_4_vect = (1:number_of_frames)'.*dv_4_step;

type = 1;
plot_temp = plot3(rr_semiVect(end,1),rr_semiVect(end,2),rr_semiVect(end,3),'or');
for k = 2:number_of_frames
    j = j+1;
    vv_2_tras_post = vv_2_tras+dv_4_vect(k,:);
    [a_tras2, e_tras2, i_tras2, OM_tras2, om_tras2, ~] = rv2parorb(rr_2_tras, vv_2_tras_post, mu);
    [T_size_new] = sizeOfTime(T_size,orbit_duration,a_tras2);
    [theta_vect] = calculateThetaVect(mu, a_tras2, e_tras2, T_size_new);
    [rr_peri_tras2, ~] = parorb2rv(a_tras2, e_tras2, i_tras2, OM_tras2, om_tras2, theta_vect, mu);
    exportPng(rr_peri_tras2,j,filename,export,type,s,frame_rotation_angle);
end
delete(plot_temp)
scatter3(rr_semiVect(end,1),rr_semiVect(end,2),rr_semiVect(end,3), 'blue', 'filled'); % Plot maneuver point



% Plot along final orbit to reach ending point
angle1 = angle2;
angle2 = theta_f;
dt_5 = timeCalc(a_f, e_f, mu, angle1, angle2);

[T_size_new] = sizeOfTime(T_size,orbit_duration,a_f);
[theta_vect] = calculateThetaVect(mu, a_f, e_f, T_size_new);

[rr_vect_complete, vv_vect_complete] = parorb2rv(a_f,e_f,i_f,OM_f,om_f,theta_vect,mu); % complete orbit
[rr_semiVect] = semiOrb(rr_vect_complete,angle1,angle2,a_f, e_f, i_f, OM_f, om_2, mu);

type = 0;
plot_temp = plot3(rr_vect_complete(:,1),rr_vect_complete(:,2),rr_vect_complete(:,3),'--m');
for k = 1:size(rr_semiVect,1)
    j = j+1;
    rr_orb = rr_semiVect(1:k,:);
    exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle,"blu");
%     exportPng(rr_orb,j,filename,export,type,s,frame_rotation_angle,"#ff5500");
end
delete(plot_temp)


% Update results
dv_tot = dv_tot+norm(dv_4);
dt_tot = dt_tot+dt_5;



if export
    warning("Done!")
end

% Display results
disp("Delta V tot = "+dv_tot+" km/s")
fprintf("\t--> delta V1 = "+norm(dv_1)+" km/s\n")
fprintf("\t--> delta V2 = "+norm(dv_2)+" km/s\n")
fprintf("\t--> delta V3 = "+norm(dv_3)+" km/s\n")
fprintf("\t--> delta V4 = "+norm(dv_4)+" km/s\n\n")

disp("Delta t tot = "+dt_tot+" s");
fprintf("\t--> delta t1 = "+dt_1+" s\n")
fprintf("\t--> delta t2 = "+dt_2+" s\n")
fprintf("\t--> delta t3 = "+dt_3+" s\n")
fprintf("\t--> delta t4 = "+dt_4+" s\n")
fprintf("\t--> delta t5 = "+dt_5+" s\n")

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

% -------------------------------------