function [a, e, i, OM, om, theta] = rv2parorb(rr, vv, mu)

% Transformation from Cartesian state to orbital elements
%
% [a, e, i, OM, om, theta] = rv2parorb (rr, vv, mu)
%
% Input arguments:
% ----------------------------------------------------------------
% rr            [3x1]   position vector                 [km]
% vv            [3x1]   velocity vector                 [km/s]
% mu            [1x1]   gravitational parameter         [km^3/s^2]
% 
% Output arguments:
% -----------------------------------------------------------------
% a             [1x1]   semi-major axis                 [km]
% e             [1x1]   eccentricity                    [-]
% i             [1x1]   inclination                     [rad]
% OM            [1x1]   RAAN                            [rad]
% om            [1x1]   pericenter anomaly              [rad]
% theta         [1x1]   true anomaly                    [rad]

I = [1 0 0];
J = [0 1 0];
K = [0 0 1];

r = norm(rr);                                       % Scalar distance 
v = norm(vv);                                       % Scalar Velocity

vr = dot(rr, vv)/r;                                 % Radial velocity

E = v^2/2 - mu/r;                                   % Specific mechanical energy
a = -mu/(2*E);                                      % Semi-major axis

h = cross(rr,vv);                                   % Specific angular momentum

i = Angle(acos(h(3)/norm(h)));                      % Inclination

N = cross(K, h);                                    % Nodes line

if (N(2) < 0)
    OM = Angle((2*pi) - acos(N(1)/norm(N)));
else
    OM = Angle(acos(N(1)/norm(N)));
end

ee = cross(vv, h)/mu - rr/r;
e = norm(ee);

if (ee(3) < 0)
    om = Angle(2*pi - acos(dot(N,ee)/(norm(N)*e)));
else
    om = Angle(acos(dot(N,ee)/(norm(N)*e)));
end

if (vr < 0)
    theta = Angle(2*pi - real(acos(dot(ee,rr)/(r*e))));
else
    theta = Angle(acos(dot(ee,rr)/(r*e)));
end

end