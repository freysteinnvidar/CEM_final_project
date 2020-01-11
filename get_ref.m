function [n,pos] = get_ref(WG0,l_rwg,R,theta,rho,RL,face)
global dx
L_tot = WG0(1) + l_rwg + R;
R = R*dx;
rho = rho*dx;
Diff = abs(rho - R/cos(theta));
RR = dx*RL/cos(theta);

A = [RR^3 RR^2;RR^3/3 RR^2/2];
B = [0;Diff];
X = linsolve(A,B);
n = 1 + X(1).*(linspace(0,RR,80)).^2 + X(2).*(linspace(0,RR,80));
pos = [dx.*(linspace(L_tot - 41,L_tot,80));face*dx +  dx*linspace(1,RL,80)*tan(theta)];
