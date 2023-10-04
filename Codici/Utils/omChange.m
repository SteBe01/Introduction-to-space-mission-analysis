function [dv, theta_3_a, theta_3_b, om_f] = omChange(a, e, om_i, om_f, mu)

% Transform the orbit changing the pericenter (om e theta)
%
% [dv, theta_3_a, theta_3_b, w_f] = omChange(a, e, om_i, om_f, mu)
%
% Input arguments:
% ----------------------------------------------------------------
% a             [1x1]   semi-major axis                 [Km]
% e             [1x1]   eccentricity                    [-]
% om_i          [1x1]   pericenter anomaly (initial)    [rad]
% om_f          [1x1]   pericenter anomaly (final)      [rad]
% mu            [1x1]   gravitational parameters        [Km^3/s^2]
% 
% -----------------------------------------------------------------
% Output arguments:
% 
% dv            [1x1]   delta v                         [km/s]
% theta_3_a     [1x1]   encounter point 1               [rad]
% theta_3_b     [1x1]   encounter point 2               [rad]
% om_f          [1x1]   pericenter anomaly (updated)    [rad]

    dw=om_f-om_i;

    if dw==0
        error("dw == 0")
    end

    if dw>pi/2 && dw<(3*pi)/2
        dw=dw-pi;
        om_f = om_f-pi;
    end
    
    theta_2_a=dw/2;
    theta_2_b=dw/2+pi;

    theta_3_a=2*pi-theta_2_a;
    theta_3_b=2*pi-theta_2_b;

    p = a*(1-e^2);
    dv = 2*sqrt(mu/p)*e*sin(dw/2);

end