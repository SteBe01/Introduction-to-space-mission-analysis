function [rr_vect_new] = semiOrb(rr_vect_complete,theta1,theta2,a,e,i,OM,om,mu)

% Cuts and extends (if asked) the complete orbit
%
% [rr_vect_new] = semiOrb(rr_vect_complete,theta1,theta2,a, e, i, OM, om, mu)
%
% Input arguments:
% ----------------------------------------------------------------
% rr_vect_complete  [nx3]   complete vector                 [km]
% theta1            [1x1]   starting point                  [rad]
% theta2            [1x1]   end point                       [rad]
% a                 [1x1]   semi-major axis                 [km]
% e                 [1x1]   eccentricity                    [-]
% i                 [1x1]   inclination                     [rad]
% OM                [1x1]   RAAN                            [rad]
% om                [1x1]   pericenter anomaly              [rad]
% mu                [1x1]   gravitational parameters        [Km^3/s^2]
% 
% Output arguments:
% -----------------------------------------------------------------
% rr_vect_new       [nx3]   final vector                    [km]

if theta1==0 && theta2==0
    warning("theta1 = theta2 = 0!");
end

roundoffproblem=0;
if theta2>=theta1
    while theta1>=2*pi
        theta1=theta1-2*pi;
        theta2=theta2-2*pi;
        roundoffproblem=1;
    end
end

if theta1>theta2
    [vect1] = semiOrbInside(rr_vect_complete,theta1,Angle(2*pi));
    vect1=vect1((1:end-1),:);
    [vect2] = semiOrbInside(rr_vect_complete,Angle(0),theta2);
    rr_vect_new=[vect1; vect2];
elseif theta2-theta1<2*pi
    if theta2<2*pi
        [rr_vect_new] = semiOrbInside(rr_vect_complete,theta1,theta2);
    else
        [vect1] = semiOrbInside(rr_vect_complete,theta1,Angle(2*pi));
        vect1=vect1((1:end-1),:);
        [vect2] = semiOrbInside(rr_vect_complete,Angle(0),theta2-2*pi);
        rr_vect_new=[vect1; vect2];
    end
else
    if theta1>0
        [vect1] = semiOrbInside(rr_vect_complete,theta1,Angle(2*pi));
        vect1=vect1((1:end-1),:);
        dummy=vect1;
        theta2=theta2-2*pi;
    else
        dummy=[];
    end
    iter=floor((theta2)/(2*pi));
    if iter
        for j=1:iter
            [vect] = semiOrbInside(rr_vect_complete,Angle(0),Angle(2*pi));
            vect=vect((1:end-1),:);
            dummy=[dummy; vect];
        end
        theta2=theta2-2*pi*iter;
    end
    [vect2] = semiOrbInside(rr_vect_complete,Angle(0),theta2);
    rr_vect_new=[dummy; vect2];
end

function [vect] = semiOrbInside(rr_vect_complete,theta1,theta2)

    [radius1, ~] = parorb2rv(a, e, i, OM, om, theta1, mu);
    [radius2, ~] = parorb2rv(a, e, i, OM, om, theta2, mu);
    
    stop1=1;
    i=0;
    dim=length(rr_vect_complete);
    if theta1>=pi
        i=floor(dim/2);
    end
    
    while stop1==1
        if theta2-theta1<=pi
            i=i+1;
            if theta1<pi
                if norm(rr_vect_complete(i,:))>=norm(radius1)
                    stop1=0;
                end
            else
                if norm(rr_vect_complete(i,:))<=norm(radius1)
                    stop1=0;
                end
            end
        else
            i=i+1;
            if norm(rr_vect_complete(i,:))>=norm(radius1)
                stop1=0;
            end
        end
        if i>=dim && roundoffproblem
            stop1=0;
            warning("Result may be affected by round off problems");
        end
        if i>= size(rr_vect_complete,1)
            stop1=0;
        end
    end
    
    
    if theta2<pi
        stop=1;
        k=1;
        while stop
            vect(k,:)=rr_vect_complete(k+i-1,:);
            if norm(vect(k,:))>norm(radius2)
                stop=0;
            end
            k=k+1;
            if k==ceil(dim/2)+1
                break
            end
        end
    elseif round(theta2.value-pi,5)>0 && round(theta2.value-2*pi,5)<0
        stop1=1;
        stop2=1;
        k=1;
        while stop1==1 || stop2==1
            vect(k,:)=rr_vect_complete(k+i-1,:);
            if norm(vect(k,:))>norm(radius2)
                stop1=0;
            end
            if stop1==0 && norm(vect(k,:))<norm(radius2)
                stop2=0;
            end
            k=k+1;
            if k==dim+1-i
                break
            end
        end
    elseif theta2==pi
        vect=rr_vect_complete((i:ceil(dim/2)),:);
    else
        vect=rr_vect_complete((i:dim),:);
    end

end


end