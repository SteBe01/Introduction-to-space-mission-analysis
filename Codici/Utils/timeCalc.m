function dt_12 = timeCalc(a, e, mu, theta_1, theta_2)

% Calculates the orbital time between two points
%
% dt_12 = timeCalc(a, e, mu, theta_1, theta_2)
%
% Input arguments:
% ----------------------------------------------------------------
% a             [1x1]   semi-major axis                 [Km]
% e             [1x1]   eccentricity                    [-]
% mu            [1x1]   gravitational parameters        [Km^3/s^2]
% theta_1       [1x1]   first angle                     [rad]
% theta_2       [1x1]   second angle                    [rad]
% 
% -----------------------------------------------------------------
% Output arguments:
% 
% dt            [1x1]   delta t                         [s]

    th_mult = sqrt((1-e)/(1+e));
    E1 = Angle(2*atan(th_mult*tan(theta_1/2)));
    E2 = Angle(2*atan(th_mult*tan(theta_2/2)));

    t_mult = sqrt(a^3/mu);
    t1 = t_mult*(double(E1) - e*sin(E1));
    t2 = t_mult*(double(E2) - e*sin(E2));

    if theta_1 > theta_2
        T_orb = 2*pi*sqrt(a^3/mu);
        dt_12 = (T_orb - t1) + t2;
    else
        dt_12 = t2 - t1;
    end
    
end