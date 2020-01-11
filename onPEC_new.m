function [isPEC,b1,b2,b21,b22] = onPEC_new(ie,je,a,b,l_rwg,R,WG0,horn,wg,outer)
isPEC = zeros(ie,je);
x0 = WG0(1);
y0 = WG0(2);
t = 6;
% PEC of the waveguide
isPEC(x0:(x0+l_rwg),[ceil(y0-a/2)-t/2:ceil(y0 - a/2), ceil(y0 + a/2):ceil(y0 + a/2)+t/2]) = 1;
if wg
    isPEC(x0-t:x0,ceil(y0-a/2)-t/2:ceil(y0 + a/2)+t/2) = 1;
end
%PEC of the horn
if horn
    slope = 2*R/(b-a);
    slope_cell = floor(slope);
    mod_slope = ceil(1/(slope - slope_cell));
    horn_max = R*slope;
    pos_horn = x0+l_rwg;
    count = 1;
    fw = round(2*R/slope_cell) + a - b;
    
    for horn = 1:R
        if mod(horn,slope_cell) == 0 && mod(horn,2*round(R/fw)) ~= 0
            count = count + 1;
        end
        isPEC(pos_horn + horn,y0 + a/2 + count:y0 + a/2 + count + t/2) = 1;
        isPEC(pos_horn + horn,y0 - a/2 - count - t/2:y0 - a/2 - count) = 1;
        if horn == R-2*20
            b22 = y0 + a/2 + count;
            b21 = y0 - a/2 - count;
        end
    end
end
b2 =   y0 + a/2 + count;
b1 = y0 - a/2 - count;
% small flare
% isPEC(pos_horn + horn:pos_horn + horn+2*t,y0 + a/2 + count:y0 + a/2 + count + t/2) = 1;
% isPEC(pos_horn + horn:pos_horn + horn+2*t,y0 - a/2 - count - t/2:y0 - a/2 - count) = 1;
% Adding PEC to the outermost layer
if outer ~= 0
    isPEC(1:outer,:) = 1;
    isPEC(end-outer-1:end,:) = 1;
    isPEC(:,1:outer) = 1;
    isPEC(:,end-outer-1:end) = 1;
end

