function [a, e, i, OM, om, theta1, theta2] = twoPointsOrbit(p1,p2,choice)

% Calculates the orbit between two given points (mu is internal)
% grater(p1/p2) = apocenter point
%
% [a, e, i, OM, om, theta1, theta2] = twoPointsOrbit(p1,p2,choice)
%
% Input arguments:
% ----------------------------------------------------------------
% p1            [1x3]   first position                          [km]
% p2            [1x3]   second position                         [km]
% choice        [1x1]   direction of the plane (1=inverted)     [-]
% 
% -----------------------------------------------------------------
% Output arguments:
% 
% a             [1x1]   semi-major axis                         [km]
% e             [1x1]   eccentricity                            [-]
% i             [1x1]   inclination                             [rad]
% OM            [1x1]   RAAN                                    [rad]
% om            [1x1]   pericenter anomaly                      [rad]
% theta1        [1x1]   angle of the first point                [rad]
% theta2        [1x1]   angle of the second point               [rad]

invert=0;
if norm(p2)>norm(p1)
    dummy=p1;
    p1=p2;
    p2=dummy;
    invert=1;
end

p1_norm=p1/norm(p1);
p2_norm=p2/norm(p2);

h=cross(p1_norm,p2_norm);
h_norm=h/norm(h);

if choice
    h_norm=-h_norm; % flips the orbit
end

n=cross([0 0 1],h_norm);
n_norm=n/norm(n);

if dot(n_norm,[0 1 0])>=0
    OM=Angle(acos(dot([1 0 0],n_norm)));
else
    OM=2*pi-Angle(acos(dot([1 0 0],n_norm)));
end

e_norm=-p1_norm;

if dot(e_norm,[0 0 1])>=0
    om=Angle(acos(dot(n_norm,e_norm)));
else
    om=2*pi-Angle(acos(dot(n_norm,e_norm)));
end

i=Angle(acos(dot(h_norm,[0 0 1])));


Rot_OM = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];    % Rotation around z axis with angle OM
Rot_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];         % Rotation around x axis with angle i
Rot_om = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];    % Rotation around z axis with angle om

Rot_matrix = Rot_om* Rot_i*Rot_OM;
p1_rot=Rot_matrix*p1';
p2_rot=Rot_matrix*p2';

theta1=Angle(atan2(p1_rot(1),p1_rot(2)))+3*pi/2;
theta2=Angle(atan2(p2_rot(1),p2_rot(2)))+3*pi/2;

e=(norm(p2_rot)-norm(p1_rot))/(norm(p1_rot)*cos(theta1)-norm(p2_rot)*cos(theta2));
a=(norm(p1_rot)+norm(p1_rot)*e*cos(theta1))/(1-e^2);

mu=398600;
[R2, ~] = parorb2rv(a, e, i, OM, om, theta2, mu);
if round(norm(R2-p2),5)~=0
    theta2=-theta2;     % handled as "Angle" class
end

if invert
    dummy=theta1;
    theta1=theta2;
    theta2=dummy;
end

end