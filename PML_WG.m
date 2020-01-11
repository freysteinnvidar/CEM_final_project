function [Sigmax,isPEC] = PML_WG(Sigmax,WG0,n,a,PML,S_WG,isPEC,opt)
x0 = WG0(1);
y0 = WG0(2);
S_x = S_WG.*((1:PML)./PML).^n;
count = 0;
t = 0;
tt = 3;
depth = 19;
if opt
    for x = x0-PML+depth+1:x0+depth
        Sigmax(x,WG0(2)-a/2-t:WG0(2)+a/2+t) = S_x(end-count).*ones(1,a+1+2*t);
        count = count + 1;
    end
else
    isPEC(x0-tt:x0,y0-a/2-tt:y0 + a/2+tt) =1;
end

