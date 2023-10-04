function [dv_1, dv_2, theta_1, theta_2, om_2] = planeChange(a, e, i_1, OM_1, om_1, i_2, OM_2, mu)
    
% Calculates the plane change
%
% [dv_1, dv_2, theta_1, theta_2, om_2] = planeChange(a, e, i_1, OM_1, om_1, i_2, OM_2, mu)
%
% Input arguments:
% ----------------------------------------------------------------
% a             [1x1]   semi-major axis                 [km]
% e             [1x1]   eccentricity                    [-]
% i_1           [1x1]   inclination                     [rad]
% OM_1          [1x1]   RAAN                            [rad]
% om_1          [1x1]   pericenter anomaly              [rad]
% i_2           [1x1]   inclination                     [rad]
% OM_2          [1x1]   RAAN                            [rad]
% mu            [1x1]   gravitational parameters        [km^3/s^2]
% 
% Output arguments:
% -----------------------------------------------------------------
% dv_1          [1x3]   delta v 1                       [km/s]
% dv_2          [1x3]   delta v 2                       [km/s]
% theta_1       [1x1]   first angle of plane change     [rad]
% theta_2       [1x1]   second angle of plane change    [rad]

    % delta OMEGA and delta i
    d_OM = OM_2-OM_1;
    d_i = i_2 - i_1;
    alpha = Angle(acos(cos(i_1)*cos(i_2)+sin(i_1)*sin(i_2)*cos(d_OM)));
    
    % cos(u1) and cos(u2)
    c_u1 = -(cos(i_2) - cos(alpha)*cos(i_1))/(sin(alpha)*sin(i_1));
    c_u2 = (cos(i_1) - cos(alpha)*cos(i_2))/(sin(alpha)*sin(i_2));

    % sin(u1) and sin(u2)
    s_u1 = sin(d_OM)/sin(alpha)*sin(i_2);
    s_u2 = sin(d_OM)/sin(alpha)*sin(i_1);

    u1 = atan2(s_u1, c_u1);
    u2 = atan2(s_u2, c_u2);

    if ~d_OM.wasNegative && ~d_i.wasNegative || d_OM.wasNegative && d_i.wasNegative
        theta_1 = u1 - om_1;
        theta_2 = theta_1;
        om_2 = u2 - theta_2;
    elseif d_OM.wasNegative && ~d_i.wasNegative|| ~d_OM.wasNegative && d_i.wasNegative
        theta_1 = u1 - om_1 + pi;
        theta_2 = theta_1;
        om_2 = u2 - theta_2 + pi;
    end

    theta_2 = theta_1+pi;

    p = a*(1-e^2);
    v_theta_1 = sqrt(mu/p)*(1+e*cos(theta_1));
    dv_1 = abs(2*v_theta_1*sin(alpha/2));

    v_theta_f_2 = sqrt(mu/p)*(1+e*cos(theta_2));
    dv_2 = abs(2*v_theta_f_2*sin(alpha/2));

end