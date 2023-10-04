function [dv_tot, dt_tot, dv_vect, dt_vect] = biellipticCost(R_p_tras_1, V_i, r_a_m, R_p_tras_2, V_f, i_i, OM_i, om_i, mu, invert, choice, orb_plot)
% This function was created to help during the biellipticTest.m development.
% The return values are total cost and time and the vectors containing
% cost and time for each maneuver

    dv_vect = []; dt_vect = [];

    r_p_tras_1 = norm(R_p_tras_1); 

    % First orbit (same plane as initial's one)
    a_tras_1 = (r_p_tras_1 + r_a_m)/2;
    e_tras_1 = (r_a_m - r_p_tras_1)/(r_a_m + r_p_tras_1);
    if invert
        om_tras_1 = om_i + pi;
    else
        om_tras_1 = om_i;
    end

    % Change from initial to elliptic
    [~, V_p_tras_1] = parorb2rv(a_tras_1, e_tras_1, i_i, OM_i, om_i, Angle(0), mu);
    dv_1 = norm(V_p_tras_1 - V_i);
    dv_vect = [dv_vect; dv_1];

    % Time from 1st transfer's perigee to apogee (where second maneuver
    % happens)
    dt_p_a_1 = timeCalc(a_tras_1, e_tras_1, mu, Angle(0), Angle(pi));
    dt_vect = [dt_vect; dt_p_a_1];

    % Plane change
    [R_a_m_1, V_a_m_1] = parorb2rv(a_tras_1, e_tras_1, i_i, OM_i, om_tras_1, Angle(pi), mu);
    
    % Second orbit
    [a_tras_2, e_tras_2, i_tras_2, OM_tras_2, om_tras_2, th_1, th_2] = twoPointsOrbit(R_a_m_1, R_p_tras_2, choice);
    [R_a_m_2, V_a_m_2] = parorb2rv(a_tras_2, e_tras_2, i_tras_2, OM_tras_2, om_tras_2, th_1, mu);

    if abs(norm(R_a_m_1 - R_a_m_2)) > 1e-6 || abs(norm(R_a_m_2) - r_a_m) > 1e-6
        error("The calculated radius do not coincide!");
    end

    dv_2 = norm(V_a_m_2 - V_a_m_1);
    dv_vect = [dv_vect; dv_2];

    % Change from elliptic to final
    [~, V_p_tras_2] = parorb2rv(a_tras_2, e_tras_2, i_tras_2, OM_tras_2, om_tras_2, th_2, mu);
    dv_3 = norm(V_f -  V_p_tras_2);
    dv_vect = [dv_vect; dv_3];

    % Time from 2nd transfer's apogee to perigee
    dt_a_p_2 = timeCalc(a_tras_2, e_tras_2, mu, th_1, th_2);
    dt_vect = [dt_vect; dt_a_p_2];

    % Total time and velocity variation needed
    dv_tot = sum(abs(dv_vect));
    dt_tot = sum(dt_vect);

    % Plot
    if orb_plot
        plotOrbit(a_tras_1, e_tras_1, i_i, OM_i, om_tras_1, 0, pi, 0.01, mu, 'c');
        plotOrbit(a_tras_2, e_tras_2, i_tras_2, OM_tras_2, om_tras_2, th_1, th_2, 0.01, mu, 'm');
    end

end