%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is used different configuration of the bielliptic maneuver.  %
% The results from this file are then as the data needed to generate the %
% animation with biellipticAnimation.m                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clear; close all; clc;
addpath("..\..\Plot\","..\..\Utils\");

% initial orbit
R_i = [-1788.3462 -9922.9190 -1645.8335];
V_i = [5.6510 -1.1520 -1.8710];
mu = 398600;
[a_i, e_i, i_i, OM_i, om_i, theta_i] = rv2parorb(R_i, V_i, mu);

% final orbit
a_f = 13290;
e_f = 0.3855;
i_f = Angle(0.9526);
OM_f = Angle(2.5510);
om_f = Angle(2.2540);
theta_f = Angle(3.0360);
[R_f, V_f] = parorb2rv(a_f, e_f, i_f, OM_f, om_f, theta_f, mu);

% Apogee radii
r_a_m_vect = (20:90)*1e3;

% Results data vector
dv_test = [];
dt_test = [];
dv_vect_test = [];
dt_vect_test = [];

orb_plot = 0; % Plots all transfer orbits.
% Note: This plot is meant as a debugging feature, it's not advised to use
% it, unless needed
if orb_plot
    figure("WindowState","maximized");
end
%% Test 1: Initial Perigee - Final Perigee

% Plot
if orb_plot
    subplot(3, 3, 1);
    plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, 0.01, mu, 'b');
    hold on; grid on; axis equal;
    plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, 0.01, mu, 'r');
    earthPlot;
    title("Test 1");
end

% Initial orbit's perigee chosen as bielliptic starting point
[R_p_i, V_p_i] = parorb2rv(a_i, e_i, i_i, OM_i, om_i, Angle(0), mu);

% Final orbit's perigee chosen as bielliptic ending point
[R_p_f, V_p_f] = parorb2rv(a_f, e_f, i_f, OM_f, om_f, Angle(0), mu);

dv_tot_vect = [];
dt_tot_vect = [];

for r_a_m = r_a_m_vect
    [dv_tot, dt_tot] = biellipticCost(R_p_i, V_p_i, r_a_m, R_p_f, V_p_f, i_i, OM_i, om_i, mu, 0, 1, orb_plot);
    dv_tot_vect = [dv_tot_vect; dv_tot];
    dt_tot_vect = [dt_tot_vect; dt_tot];
end

dt_tot_vect = dt_tot_vect + timeCalc(a_i, e_i, mu, theta_i, Angle(0)) + ...
                timeCalc(a_f, e_f, mu, Angle(0), theta_f);

dv_test = [dv_test dv_tot_vect];
dt_test = [dt_test dt_tot_vect];

%% Test 2: Initial Apogee - Final Perigee

% Plot
if orb_plot
    subplot(3, 3, 2);
    plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, 0.01, mu, 'b');
    hold on; grid on; axis equal;
    plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, 0.01, mu, 'r');
    earthPlot;
    title("Test 2");
end

% Initial orbit's apogee chosen as bielliptic starting point
[R_a_i, V_a_i] = parorb2rv(a_i, e_i, i_i, OM_i, om_i, Angle(pi), mu);

% Final orbit's perigee chosen as bielliptic ending point
[R_p_f, V_p_f] = parorb2rv(a_f, e_f, i_f, OM_f, om_f, Angle(0), mu);

dv_tot_vect = [];
dt_tot_vect = [];

for r_a_m = r_a_m_vect
    [dv_tot, dt_tot] = biellipticCost(R_a_i, V_a_i, r_a_m, R_p_f, V_p_f, i_i, OM_i, om_i, mu, 1, 1, orb_plot);
    dv_tot_vect = [dv_tot_vect; dv_tot];
    dt_tot_vect = [dt_tot_vect; dt_tot];
end
dt_tot_vect = dt_tot_vect + timeCalc(a_i, e_i, mu, theta_i, Angle(pi)) + ...
                timeCalc(a_f, e_f, mu, Angle(0), theta_f);

dv_test = [dv_test dv_tot_vect];
dt_test = [dt_test dt_tot_vect];

%% Test 3: Starting Point - Final perigee

% Plot
if orb_plot
    subplot(3, 3, 3);
    plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, 0.01, mu, 'b');
    hold on; grid on; axis equal;
    plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, 0.01, mu, 'r');
    earthPlot;
    title("Test 3");
end

% Starting point chosen as bielliptic starting point

% Final orbit's perigee chosen as bielliptic ending point
[R_p_f, V_p_f] = parorb2rv(a_f, e_f, i_f, OM_f, om_f, Angle(0), mu);

dv_tot_vect = [];
dt_tot_vect = [];

for r_a_m = r_a_m_vect
    [dv_tot, dt_tot] = biellipticCost(R_i, V_i, r_a_m, R_p_f, V_p_f, i_i, OM_i, om_i, mu, 1, 1, orb_plot);
    dv_tot_vect = [dv_tot_vect; dv_tot];
    dt_tot_vect = [dt_tot_vect; dt_tot];
end

dt_tot_vect = dt_tot_vect + timeCalc(a_f, e_f, mu, Angle(0), theta_f);

dv_test = [dv_test dv_tot_vect];
dt_test = [dt_test dt_tot_vect];

%% Test 4: Initial Perigee - Final Apogee

% Plot
if orb_plot
    subplot(3, 3, 4);
    plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, 0.01, mu, 'b');
    hold on; grid on; axis equal;
    plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, 0.01, mu, 'r');
    earthPlot;
    title("Test 4");
end

% Initial orbit's perigee chosen as bielliptic starting point
[R_p_i, V_p_i] = parorb2rv(a_i, e_i, i_i, OM_i, om_i, Angle(0), mu);

% Final orbit's apogee chosen as bielliptic ending point
[R_a_f, V_a_f] = parorb2rv(a_f, e_f, i_f, OM_f, om_f, Angle(pi), mu);

dv_tot_vect = [];
dt_tot_vect = [];

for r_a_m = r_a_m_vect
    [dv_tot, dt_tot] = biellipticCost(R_p_i, V_p_i, r_a_m, R_a_f, V_a_f, i_i, OM_i, om_i, mu, 0, 0, orb_plot);
    dv_tot_vect = [dv_tot_vect; dv_tot];
    dt_tot_vect = [dt_tot_vect; dt_tot];
end

dt_tot_vect = dt_tot_vect + timeCalc(a_i, e_i, mu, theta_i, Angle(0)) + ...
                timeCalc(a_f, e_f, mu, Angle(pi), theta_f);

dv_test = [dv_test dv_tot_vect];
dt_test = [dt_test dt_tot_vect];

%% Test 5: Initial Apogee - Final Apogee

% Plot
if orb_plot
    subplot(3, 3, 5);
    plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, 0.01, mu, 'b');
    hold on; grid on; axis equal;
    plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, 0.01, mu, 'r');
    earthPlot;
    title("Test 5");
end

% Initial orbit's apogee chosen as bielliptic starting point
[R_a_i, V_a_i] = parorb2rv(a_i, e_i, i_i, OM_i, om_i, Angle(pi), mu);

% Final orbit's apogee chosen as bielliptic ending point
[R_a_f, V_a_f] = parorb2rv(a_f, e_f, i_f, OM_f, om_f, Angle(pi), mu);

dv_tot_vect = [];
dt_tot_vect = [];

for r_a_m = r_a_m_vect
    [dv_tot, dt_tot] = biellipticCost(R_a_i, V_a_i, r_a_m, R_a_f, V_a_f, i_i, OM_i, om_i, mu, 1, 0, orb_plot);
    dv_tot_vect = [dv_tot_vect; dv_tot];
    dt_tot_vect = [dt_tot_vect; dt_tot];
end

dt_tot_vect = dt_tot_vect + timeCalc(a_i, e_i, mu, theta_i, Angle(pi)) + ...
                timeCalc(a_f, e_f, mu, Angle(pi), theta_f);

dv_test = [dv_test dv_tot_vect];
dt_test = [dt_test dt_tot_vect];

%% Test 6: Starting Point - Final Apogee

% Plot
if orb_plot
    subplot(3, 3, 6);
    plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, 0.01, mu, 'b');
    hold on; grid on; axis equal;
    plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, 0.01, mu, 'r');
    earthPlot;
    title("Test 6");
end

% Starting point chosen as bielliptic starting point

% Final orbit's apogee chosen as bielliptic ending point
[R_a_f, V_a_f] = parorb2rv(a_f, e_f, i_f, OM_f, om_f, Angle(pi), mu);

dv_tot_vect = [];
dt_tot_vect = [];

for r_a_m = r_a_m_vect
    [dv_tot, dt_tot] = biellipticCost(R_i, V_i, r_a_m, R_a_f, V_a_f, i_i, OM_i, om_i, mu, 1, 0, orb_plot);
    dv_tot_vect = [dv_tot_vect; dv_tot];
    dt_tot_vect = [dt_tot_vect; dt_tot];
end

dt_tot_vect = dt_tot_vect + timeCalc(a_f, e_f, mu, Angle(pi), theta_f);

dv_test = [dv_test dv_tot_vect];
dt_test = [dt_test dt_tot_vect];

%% Test 7: Initial Perigee - Ending Point

% Plot
if orb_plot
    subplot(3, 3, 7);
    plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, 0.01, mu, 'b');
    hold on; grid on; axis equal;
    plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, 0.01, mu, 'r');
    earthPlot;
    title("Test 7");
end

% Initial orbit's perigee chosen as bielliptic starting point
[R_p_i, V_p_i] = parorb2rv(a_i, e_i, i_i, OM_i, om_i, Angle(0), mu);

% Ending point chosen as bielliptic ending point

dv_tot_vect = [];
dt_tot_vect = [];
dv_vect_vect = [];
dt_vect_vect = [];

for r_a_m = r_a_m_vect
    [dv_tot, dt_tot, dv_vect, dt_vect] = biellipticCost(R_p_i, V_p_i, r_a_m, R_f, V_f, i_i, OM_i, om_i, mu, 0, 0, orb_plot);
    dv_tot_vect = [dv_tot_vect; dv_tot];
    dt_tot_vect = [dt_tot_vect; dt_tot];
    dv_vect_vect = [dv_vect_vect dv_vect];
    dt_vect_vect = [dt_vect_vect dt_vect];
end

dt_tot_vect = dt_tot_vect + timeCalc(a_i, e_i, mu, theta_i, Angle(0));

dv_test = [dv_test dv_tot_vect];
dt_test = [dt_test dt_tot_vect];
dv_vect_test = [dv_vect_test dv_vect_vect];
dt_vect_test = [dt_vect_test dt_vect_vect];

%% Test 8: Initial Apogee - Ending Point

% Plot
if orb_plot
    subplot(3, 3, 8);
    plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, 0.01, mu, 'b');
    hold on; grid on; axis equal;
    plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, 0.01, mu, 'r');
    earthPlot;
    title("Test 8");
end

% Initial orbit's apogee chosen as bielliptic starting point
[R_a_i, V_a_i] = parorb2rv(a_i, e_i, i_i, OM_i, om_i, Angle(pi), mu);

% Ending point chosen as bielliptic ending point

dv_tot_vect = [];
dt_tot_vect = [];

for r_a_m = r_a_m_vect
    [dv_tot, dt_tot] = biellipticCost(R_a_i, V_a_i, r_a_m, R_f, V_f, i_i, OM_i, om_i, mu, 1, 0, orb_plot);
    dv_tot_vect = [dv_tot_vect; dv_tot];
    dt_tot_vect = [dt_tot_vect; dt_tot];
end

dt_tot_vect = dt_tot_vect + timeCalc(a_i, e_i, mu, theta_i, Angle(pi));

dv_test = [dv_test dv_tot_vect];
dt_test = [dt_test dt_tot_vect];

%% Test 9: Starting Point - Ending Point

% Plot
if orb_plot
    subplot(3, 3, 9);
    plotOrbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, 0.01, mu, 'b');
    hold on; grid on; axis equal;
    plotOrbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, 0.01, mu, 'r');
    earthPlot;
    title("Test 9");
end

% Starting point chosen as bielliptic starting point

% Ending point chosen as bielliptic ending point

dv_tot_vect = [];
dt_tot_vect = [];

for r_a_m = r_a_m_vect
    [dv_tot, dt_tot] = biellipticCost(R_i, V_i, r_a_m, R_f, V_f, i_i, OM_i, om_i, mu, 1, 0, orb_plot);
    dv_tot_vect = [dv_tot_vect; dv_tot];
    dt_tot_vect = [dt_tot_vect; dt_tot];
end

dv_test = [dv_test dv_tot_vect];
dt_test = [dt_test dt_tot_vect];

%% Plot results

figure("WindowState","maximized");
subplot(1, 2, 1);
plot(r_a_m_vect, dv_test, 'LineWidth', 1.25);
grid on;
title("Total \Deltav");
xlabel("Apogee radius");
legend("Test 1", "Test 2", "Test 3", "Test 4", "Test 5", "Test 6", "Test 7", "Test 8", "Test 9",...
        'Location', 'best');
subplot(1, 2, 2);
plot(r_a_m_vect, dt_test, 'LineWidth', 1.25);
grid on;
title("Title \Deltat for the transfer");
xlabel("Apogee radius");
legend("Test 1", "Test 2", "Test 3", "Test 4", "Test 5", "Test 6", "Test 7", "Test 8", "Test 9",...
        'Location', 'northwest');

%% Conclusion

% NOTE: From the plots we can see that the Tests 1,2,3 generate maneuvers that
% pass through Earth, so all three tests will be discarded and not further
% considered in the results analysis.
%
% If we analyze the different dv needed for the maneuvers, we can notice
% that the lowest are given by Test 4 and Test 7.
% Between the two we decided to chose Test 7 as the source of data for
% biellipticAnimation.m
% This decision was made noticing that this was the only test where the dv reached 
% a minimum and then started climbing again.
% Therefor we can calculate the radius associated with the minimum
% velocity:

test_n = 7;

[val, pos] = min(dv_test(:, test_n));
r_a_m_result = r_a_m_vect(pos)
dv_result = dv_test(pos, test_n)
dt_result = dt_test(pos, test_n)/60^2

h=figure;
plot(r_a_m_vect, dv_test(:, test_n), 'b','LineWidth', 1.5);
hold on, grid on

xlabel('Apocenter radius [Km]',FontSize=15);
ylabel('\Deltav [Km/s]',FontSize=15);

% title("Test " + test_n);
scatter(r_a_m_result, dv_result, 'green', 'filled');
plot(r_a_m_vect, dv_vect_test(2, :), '--r', 'LineWidth',1);
plot(r_a_m_vect, dv_vect_test(1,:) + dv_vect_test(3,:), '--m', 'LineWidth', 1);
legend("\Deltav", "Minimum point", "Plane change", "Orbit change",'Location','east');


% other figure

figure;
plot(r_a_m_vect, dv_test(:, test_n), 'b','LineWidth', 1.5);
hold on, grid on

xlabel('Apocenter radius [Km]',FontSize=15);
ylabel('\Deltav [Km/s]',FontSize=15);

scatter(r_a_m_result, dv_result, 'green', 'filled');
legend("\Deltav", "Minimum point",'Location','east');

%% export settings

% set(h,'PaperOrientation','landscape');
% set(h,'PaperPosition', [1 1 28 19]);
% print(gcf, 'Bielliptic4.pdf', '-dpdf', '-fillpage')