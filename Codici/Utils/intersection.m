function [p1,p2,angle1,angle2] = intersection(a_1,e_1,i_1,OM_1,om_1,a_2,e_2,i_2,OM_2,om_2,mu,toll,angle_inter,type)

% Calculates the intersection between two orbits
% type=0: angle_inter (arbitrary) --> no double intersection
% type=1: angle_inter --> relative to the first orbit
% type=2: angle_inter --> relative to the second orbit
%
% [p1,p2,angle1,angle2] = intersection(a_1,e_1,i_1,OM_1,om_1,a_2,e_2,i_2,OM_2,om_2,mu,toll,angle_inter,type)
%
% Input arguments:
% ----------------------------------------------------------------
% a_1           [1x1]   semi-major axis                 [km]
% e_1           [1x1]   eccentricity                    [-]
% i_1           [1x1]   inclination                     [rad]
% OM_1          [1x1]   RAAN                            [rad]
% om_1          [1x1]   pericenter anomaly              [rad]
% a_2           [1x1]   semi-major axis                 [km]
% e_2           [1x1]   eccentricity                    [-]
% i_2           [1x1]   inclination                     [rad]
% OM_2          [1x1]   RAAN                            [rad]
% om_2          [1x1]   pericenter anomaly              [rad]
% mu            [1x1]   gravitational parameters        [Km^3/s^2]
% toll          [1x1]   tollerance                      [-]
% angle_inter   [1x1]   intersection angle              [rad]
% type          [1x1]   type of angle_inter             [-]
% 
% -----------------------------------------------------------------
% Output arguments:
% 
% p1            [1x3]   first radius (angle1)           [km]
% p1            [1x3]   second radius (angle2)          [km]
% angle1        [1x1]   interseption on first orbit     [rad]
% angle2        [1x1]   interseption on second orbit    [rad]

if type==0
    angle1_1=Angle(0);
    angle1_2=Angle(0);
    angle2_1=Angle(2*pi);
    angle2_2=Angle(2*pi);
    tot_theta=angle2_1-angle1_1;
elseif type==1
    angle1_1=Angle(angle_inter-pi/2);
    angle1_2=Angle(0);
    angle2_1=Angle(angle_inter+pi/2);
    angle2_2=Angle(2*pi);
    tot_theta=angle2_1-angle1_1;
elseif type==2
    angle1_1=Angle(0);
    angle1_2=Angle(angle_inter-pi/2);
    angle2_1=Angle(2*pi);
    angle2_2=Angle(angle_inter+pi/2);
    tot_theta=angle2_2-angle1_2;
end

if angle1_1.wasNegative
    angle1_1=Angle(0);
end
if angle1_2.wasNegative
    angle1_2=Angle(0);
end
if angle2_1.wasNegative
    angle2_1=Angle(0);
end
if angle2_2.wasNegative
    angle2_2=Angle(0);
end

factor=100;
step=(tot_theta)/factor;

angle1=Angle(0);
angle2=Angle(0);
distance=500;           % dummy distance

found=0;
while step > toll
    if angle1_1.value-step<0
        angle1_1=Angle(0);
    end
    if angle1_2.value-step<0
        angle1_2=Angle(0);
    end
    if angle2_1.value-step<0
        angle2_1=Angle(0);
    end
    if angle2_2.value-step<0
        angle2_2=Angle(0);
    end

    theta_vect_1=linspace(angle1_1,angle2_1,factor)';
    theta_vect_2=linspace(angle1_2,angle2_2,factor)';
    
    [rr_1, ~] = parorb2rv(a_1, e_1, i_1, OM_1, om_1, theta_vect_1, mu);
    [rr_2, ~] = parorb2rv(a_2, e_2, i_2, OM_2, om_2, theta_vect_2, mu);
    
    for i=1:length(theta_vect_1.value)
        for k=1:length(theta_vect_2.value)
            distance_new=norm(rr_1(i,:)-rr_2(k,:));
            if distance_new<distance
                found=1;
                distance=distance_new;
                angle1=Angle(theta_vect_1.value(i));
                angle2=Angle(theta_vect_2.value(k));
            end
        end
    end
    
    angle1_1=angle1-step;
    angle2_1=angle1+step;
    angle1_2=angle2-step;
    angle2_2=angle2+step;
    
    tot_theta=2*step;
    step=(tot_theta)/factor;
end

angle1=Angle(real(angle1.value));
angle2=Angle(real(angle2.value));

[p1, ~] = parorb2rv(a_1, e_1, i_1, OM_1, om_1, angle1, mu);
[p2, ~] = parorb2rv(a_2, e_2, i_2, OM_2, om_2, angle2, mu);

if ~found
    error("No intersection found")
end

end