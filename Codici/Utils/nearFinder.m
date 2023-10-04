function [p1,p2,angle1,angle2] = nearFinder(a_1,e_1,i_1,OM_1,om_1,a_2,e_2,i_2,OM_2,om_2,mu,toll)

% Calculates the minimum distance between two orbits
%
% [p1,p2,angle1,angle2] = nearFinder(a_1,e_1,i_1,OM_1,om_1,a_2,e_2,i_2,OM_2,om_2,mu,toll)
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
% 
% -----------------------------------------------------------------
% Output arguments:
% 
% p1            [1x3]   first radius (angle1)           [km]
% p1            [1x3]   second radius (angle2)          [km]
% angle1        [1x1]   interseption on first orbit     [rad]
% angle2        [1x1]   interseption on second orbit    [rad]

angle1_1=Angle(0);
angle1_2=Angle(0);
angle2_1=Angle(2*pi);
angle2_2=Angle(2*pi);
tot_theta=angle2_1-angle1_1;
passo=(tot_theta)/100;

[rr_1, ~] = parorb2rv(a_1, e_1, i_1, OM_1, om_1, 0, mu);
[rr_2, ~] = parorb2rv(a_2, e_2, i_2, OM_2, om_2, 0, mu);

distance=norm(rr_1(1,:)-rr_2(1,:));
angle1=Angle(0);
angle2=Angle(0);

while passo > toll
    theta_vect_1=linspace(angle1_1,angle2_1,100)';
    theta_vect_2=linspace(angle1_2,angle2_2,100)';
    
    [rr_1, vv_1] = parorb2rv(a_1, e_1, i_1, OM_1, om_1, theta_vect_1, mu);
    [rr_2, vv_2] = parorb2rv(a_2, e_2, i_2, OM_2, om_2, theta_vect_2, mu);
    
    for i=1:length(theta_vect_1.value)
        for k=1:length(theta_vect_2.value)
            distance_new=norm(rr_1(i,:)-rr_2(k,:));
            if distance_new<distance
                distance=distance_new;
                [~, ~, ~, ~, ~, angle1] = rv2parorb(rr_1(i,:), vv_1(i,:), mu);
                [~, ~, ~, ~, ~, angle2] = rv2parorb(rr_2(k,:), vv_2(k,:), mu);
            end
        end
    end
    
    angle1_1=angle1-passo/2;
    angle2_1=angle1+passo/2;
    angle1_2=angle2-passo/2;
    angle2_2=angle2+passo/2;
    
    tot_theta=passo;
    passo=(tot_theta)/100;
end

angle1=Angle(real(angle1.value));
angle2=Angle(real(angle2.value));

[p1, ~] = parorb2rv(a_1, e_1, i_1, OM_1, om_1, angle1, mu);
[p2, ~] = parorb2rv(a_2, e_2, i_2, OM_2, om_2, angle2, mu);

end