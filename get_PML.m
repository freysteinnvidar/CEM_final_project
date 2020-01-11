function [Sigmax,Sigmay] = get_PML(PML,ie,je,S_max,n,Sigmax,Sigmay,shift)
% shift = 5;
S_x = S_max.*((1:PML)./PML).^n; % Growth of conductivity
% Do Right and Left where Sigmay == 0

for y = 1:je % in y
    Sigmax(end-PML+1-shift:end-shift,y) = S_x;
    Sigmax(1+shift:PML+shift,y) = fliplr(S_x);
end

for x = 1:ie % in x
    Sigmay(x,end-PML+1-shift:end-shift) = S_x;
    Sigmay(x,1+shift:PML+shift) = fliplr(S_x);
end


