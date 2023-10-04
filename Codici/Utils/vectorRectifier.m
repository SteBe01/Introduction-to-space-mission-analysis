function [rr_vect] = vectorRectifier(rr_vect)

% Calculates the correct vector (inverts half of it)
%
% [rr_vect] = vectorRectifier(rr_vect)
%
% Input arguments:
% ----------------------------------------------------------------
% rr_vect       [T_sizex1]/[1xT_size]   initial vector       [---]
% 
% -----------------------------------------------------------------
% Output arguments:
% 
% rr_vect       [T_sizex1]/[1xT_size]   output vector        [---]

inversion=0;
if size(rr_vect,1) < size(rr_vect,2)
    rr_vect = rr_vect';
    inversion = 1;
end

dim = size(rr_vect,1);

if rem(dim,2) == 0
    for i = 1:dim
        if i<=dim/2
            dummy(i,:) = rr_vect(dim/2+i,:);
        else
            dummy(i,:) = rr_vect(i-dim/2,:);
        end
    end
    dummy=[dummy; dummy(1,:)];
else
    for i = 1:dim
        if i <= floor(dim/2)
            dummy(i,:) = rr_vect(floor(dim/2)+i,:);
        else
            dummy(i,:) = rr_vect(i-floor(dim/2),:);
        end
    end
end

if inversion
    rr_vect = dummy';
else
    rr_vect = dummy;
end

end