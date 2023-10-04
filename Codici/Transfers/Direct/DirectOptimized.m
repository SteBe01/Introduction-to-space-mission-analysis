% This script generates the animation for the direct transfer between
% initial and final orbits, with optimization for dv or dt.
% The maneuver uses a transfer orbit that connects directly the starting
% and ending point.
% This script also calculates the fuel mass ratio needed for the whole transfer 
% using the Tsiolkovsky equation.
% Parameters:
%   - Optimization (line 49): 
%           -> 0: transfer with lowest dt
%           -> 1: transfer with lowest dv
%   - Export (line 54):
%           -> 0: Plot animation without saving frames
%           -> 1: Save all frame as .png images
%   - Maximize (line 53);
%           -> 0: Open figure window with default size
%           -> 1: Open figure window maximized
%

%% Animation engine

clear, clc
close all

addpath("..\..\Plot\","..\..\Utils\");

global j filename s frame_rotation_angle

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

optimization = 1;                   
T_size = 300;                       % Amount of points per orbit
I_sp = 348;                         % Engine specific impulse (Engine: Merlin 1D+ Vacuum)

maximize = 1;
export = 0;
filename = '.\Export\Direct';       % Filepath of saved frames (folder + frame name)
frame_duration = 0.0333;            % Number of seconds for each frame (in real time)
number_of_frames = 20;              % Number of frames for maneuver's animation

toll = 1e-10;                       % Tollerance used when calculating near points on orbits

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
view(70,10);
hold on, grid on, axis equal
xlabel('X [Km]',FontSize=15);
ylabel('Y [Km]',FontSize=15);
zlabel('Z [Km]',FontSize=15);

% earth
s=earthPlot;
% Pre rotation of the Earth
rotate(s, [0 0 1],280,[0 0 0])

one_orbit_time=frame_duration*T_size;
orbit_duration=2*pi*sqrt((a_i^3)/mu);

factor=one_orbit_time/orbit_duration;

earth_rotation_duration = 24*60*60+56*60+04; % 1 day == 24h 56m 04s
earth_rotation_duration=earth_rotation_duration*factor;
angular_vel=360/earth_rotation_duration;
frame_rotation_angle=angular_vel*frame_duration;

% Framing
[apo_i, ~] = parorb2rv(a_i, e_i, i_i, OM_i, om_i, pi, mu);
[apo_f, ~] = parorb2rv(a_f, e_f, i_f, OM_f, om_f, pi, mu);
max_frame=max(norm(apo_i),norm(apo_f));

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

plot3(rr_i(1),rr_i(2),rr_i(3),'xr',LineWidth=2);
[rr_f_point, ~] = parorb2rv(a_f, e_f, i_f, OM_f, om_f, theta_f, mu);
plot3(rr_f_point(1),rr_f_point(2),rr_f_point(3),'xr',LineWidth=2);



% Calculate starting and ending point
[rr_1, vv_1] = parorb2rv(a_i, e_i, i_i, OM_i, om_i, theta_i, mu);
[rr_2, vv_2] = parorb2rv(a_f, e_f, i_f, OM_f, om_f, theta_f, mu);


% The choice variable decides if orbit plane is inverted
% Generate first possible transfer orbit
choice=0;
[a_tras1, e_tras1, i_tras1, OM_tras1, om_tras1] = twoPointsOrbit(rr_1,rr_2,choice);
[~,~,~,angle_c0_1] = nearFinder(a_i,e_i,i_i,OM_i,om_i,a_tras1,e_tras1,i_tras1,OM_tras1,om_tras1,mu,toll);
[~,~,angle_c0_2,~] = nearFinder(a_tras1,e_tras1,i_tras1,OM_tras1,om_tras1,a_f,e_f,i_f,OM_f,om_f,mu,toll);
dt1 = timeCalc(a_tras1, e_tras1, mu, angle_c0_1, angle_c0_2);
[~, vv_tras_c0_i] = parorb2rv(a_tras1, e_tras1, i_tras1, OM_tras1, om_tras1, angle_c0_1, mu);
[~, vv_tras_c0_f] = parorb2rv(a_tras1, e_tras1, i_tras1, OM_tras1, om_tras1, angle_c0_2, mu);
dv_i_vect_c0=vv_tras_c0_i-vv_1;
dv_f_vect_c0=vv_2-vv_tras_c0_f;

% Generate second possible transfer orbit
choice=1;
[a_tras2, e_tras2, i_tras2, OM_tras2, om_tras2] = twoPointsOrbit(rr_1,rr_2,choice);
[~,~,~,angle_c1_1] = nearFinder(a_i,e_i,i_i,OM_i,om_i,a_tras2,e_tras2,i_tras2,OM_tras2,om_tras2,mu,toll);
[~,~,angle_c1_2,~] = nearFinder(a_tras2,e_tras2,i_tras2,OM_tras2,om_tras2,a_f,e_f,i_f,OM_f,om_f,mu,toll);
dt2 = timeCalc(a_tras2, e_tras2, mu, angle_c1_1, angle_c1_2);
[~, vv_tras_c1_i] = parorb2rv(a_tras2, e_tras2, i_tras2, OM_tras2, om_tras2, angle_c1_1, mu);
[~, vv_tras_c1_f] = parorb2rv(a_tras2, e_tras2, i_tras2, OM_tras2, om_tras2, angle_c1_2, mu);
dv_i_vect_c1=vv_tras_c1_i-vv_1;
dv_f_vect_c1=vv_2-vv_tras_c1_f;

% Calculate dv associated with the transfer orbits generated
dv_c0=norm(dv_i_vect_c0)+norm(dv_f_vect_c0);
dv_c1=norm(dv_i_vect_c1)+norm(dv_f_vect_c1);

% Chose orbit based on optimization parameter 
timeForced=0;       % If dv is equal in both points, optimize for dt
if optimization
    if dv_c0<dv_c1
        choice=0;
        angle1=angle_c0_1;
        angle2=angle_c0_2;
        dv_tot=dv_c0;
        dv_1=norm(dv_i_vect_c0);
        dv_2=norm(dv_f_vect_c0);
    elseif dv_c0>dv_c1
        choice=1;
        angle1=angle_c1_1;
        angle2=angle_c1_2;
        dv_tot=dv_c1;
        dv_1=norm(dv_i_vect_c1);
        dv_2=norm(dv_f_vect_c1);
    else
        timeForced=1;
    end
end
if ~optimization || timeForced
    if dt1<dt2
        choice=0;
        angle1=angle_c0_1;
        angle2=angle_c0_2;
        dv_tot=dv_c0;
        dv_1=norm(dv_i_vect_c0);
        dv_2=norm(dv_f_vect_c0);
    else
        choice=1;
        angle1=angle_c1_1;
        angle2=angle_c1_2;
        dv_tot=dv_c1;
        dv_1=norm(dv_i_vect_c1);
        dv_2=norm(dv_f_vect_c1);
    end
end


% Generate direct transfer orbit data based on optimization
[a_tras, e_tras, i_tras, OM_tras, om_tras] = twoPointsOrbit(rr_1,rr_2,choice);


% Animate orbit change
[i_vect,OM_vect,om_vect]=parametersSpread(number_of_frames,i_i,i_tras,OM_i,OM_tras,om_i,om_tras);
[a_vect,e_vect,~]=parametersSpread(number_of_frames,a_i,a_tras,e_i,e_tras,0,0);

type=1;
plot_temp=plot3(rr_1(1),rr_1(2),rr_1(3),'or');
for k=2:number_of_frames % without initial orbit
    j=j+1;  % global
    [T_size_new] = sizeOfTime(T_size,orbit_duration,a_vect(k));
    [theta_vect] = calculateThetaVect(mu, a_vect(k), e_vect(k), T_size_new);
    [rr_orb, ~] = parorb2rv(a_vect(k),e_vect(k),i_vect.value(k),OM_vect.value(k),om_vect.value(k),theta_vect,mu);
    exportPng(rr_orb,export,type);
end
delete(plot_temp)
scatter3(rr_1(end,1),rr_1(end,2),rr_1(end,3), 'blue', 'filled'); % Plot maneuver point
plot3(rr_i(1),rr_i(2),rr_i(3),'xr',LineWidth=2);




% Plot along initial orbit to reach maneuver point
[rr_vect_complete, ~] = parorb2rv(a_tras,e_tras,i_tras,OM_tras,om_tras,theta_vect,mu); % complete orbit
[rr_semiVect] = semiOrb(rr_vect_complete,angle1,angle2,a_tras, e_tras, i_tras, OM_tras, om_tras, mu);

type=0;
plot_temp=plot3(rr_vect_complete(:,1),rr_vect_complete(:,2),rr_vect_complete(:,3),'--m');
for k=1:size(rr_semiVect,1)
    j=j+1;
    rr_orb=rr_semiVect(1:k,:);
    exportPng(rr_orb,export,type,"blu");
%     exportPng(rr_orb,export,type,"#bf00ff");
end
delete(plot_temp)



% Animate orbit change
[i_vect,OM_vect,om_vect]=parametersSpread(number_of_frames,i_tras,i_f.value,OM_tras,OM_f,om_tras,om_f);
[a_vect,e_vect,~]=parametersSpread(number_of_frames,a_tras,a_f,e_tras,e_f,0,0);

type=1;
plot_temp=plot3(rr_semiVect(end,1),rr_semiVect(end,2),rr_semiVect(end,3),'or');
for k=2:number_of_frames % without initial orbit
    j=j+1;  % global
    [T_size_new] = sizeOfTime(T_size,orbit_duration,a_vect(k));
    [theta_vect] = calculateThetaVect(mu, a_vect(k), e_vect(k), T_size_new);
    [rr_orb, ~] = parorb2rv(a_vect(k),e_vect(k),i_vect.value(k),OM_vect.value(k),om_vect.value(k),theta_vect,mu);
    exportPng(rr_orb,export,type);
end
delete(plot_temp)
scatter3(rr_semiVect(end,1),rr_semiVect(end,2),rr_semiVect(end,3), 'blue', 'filled'); % Plot maneuver point
plot3(rr_f_point(1),rr_f_point(2),rr_f_point(3),'xr',LineWidth=2);



% Update results
[~, vv_tras_f] = parorb2rv(a_tras, e_tras, i_tras, OM_tras, om_tras, angle2, mu);
dv_2_vect=vv_2-vv_tras_f;

dt_1 = timeCalc(a_tras, e_tras, mu, angle1, angle2);
dt_tot=dt_1;



if export
    warning("Done!")
end

disp("Delta V tot = "+dv_tot+" km/s")
fprintf("\t--> delta V1 = "+norm(dv_1)+" km/s\n")
fprintf("\t--> delta V2 = "+norm(dv_2)+" km/s\n")

disp("Delta t tot = "+dt_tot+" s");
fprintf("\t--> delta t1 = "+dt_1+" s\n")

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

function exportPng(rr_vect,export,type,color)
    % Generate plot and save it as png, if export is True

    global j filename s frame_rotation_angle

    x=rr_vect(:,1);
    y=rr_vect(:,2);
    z=rr_vect(:,3);
    
    if ~type
        plot3(x,y,z,'Color',color,LineWidth=1.3);
        frame1=plot3(x(end),y(end),z(end),'or');
        rotate(s, [0 0 1],frame_rotation_angle,[0 0 0])
    else
        frame1=plot3(x,y,z,'--m');
    end

    pause(0.1)
    % Save frame as .png file
    if export
        name=append(filename,"_",num2str(j),".png");
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


function [export]=answerFunction(export)
    answer=0;
    while ~answer
        reply = input('Are you sure to create PNGs? [Y/N]:','s');
        if isempty(reply)
            reply = 'Y';
        end
        if strlength(reply)==1
            if reply=='n' || reply=='N'
                warning("'export' set to 0")
                export=0;
                answer=1;
            elseif reply=='y' || reply=='Y'
                warning("Generating PNGs...")
                answer=1;
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