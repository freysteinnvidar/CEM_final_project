%% Main section with the definition of the domain and FDTD of Maxwells equations
% Initalizing Domain and physical quantities
clear all
close all
x_cell = 400; y_cell = 260; f = 1e9; c0 = 299792458; wl = c0/f;
er0 = 8.85e-12; mur0 = 4*pi*1e-7;
global dx dt
dx = wl/20; % Physical size of each cell
dt = 0.98*dx/(c0*sqrt(2)); % Courant condition
S_air = 0; S_PEC = 6e7; % Conductivities

I_region = ones(x_cell,y_cell);
% Permiability Permitivity and magnetic conductivity
er = er0.*I_region;
mur = mur0.*I_region;
Sigmay = S_air.*I_region;
Sigmax = S_air.*I_region;

% Horn antenna dimensions With PEC boundary
% boundaries 
a =ceil(wl/2/dx); % Width of waveguide
b = ceil(wl*4/dx); % Size of aparteure
l_rwg = ceil(2*wl/dx); % length of waveuigde
R = ceil(8*wl/dx); % Length of horn antenna
WG0 = [round(4*wl/dx),round(y_cell/2)]; % Center of waveuigde
S0 = [WG0(1)+round(1/4*wl/dx),WG0(2)]; % Position of source

%Designing the dielectric lense
neff = 1.4; % Value for the hemespherical lens
flag_lens = 'none'; % hem for hemespherical lens,grad for graded index lens and anything else for no lens
shift = 0; % If one wants to terminate with PEC set shift == 2;
horn_sim = 1;
[isPEC,b1,b2,b21,b22] = onPEC_new(x_cell,y_cell,a,b,l_rwg,R,WG0,horn_sim,0,shift); % function defining pec
[er,mur] = get_lens(WG0,l_rwg,neff,er,mur,R,er0,mur0,isPEC,b1,b2,b21,b22,flag_lens); % function for getting lens
% imagesc(er) %Remove comment to see lense on grid
% Define the PML boundry
PML = 50; % Number of cells for PML
n = 7; % Order of PML - How fast it grows
S_max = 20; % Maximum value of conductivity
%PEC boundary for magnetic and electric conductivity
[Sigmax,Sigmay] = get_PML(PML,x_cell,y_cell,S_max,n,Sigmax,Sigmay,shift); % add the PML boundary
% PML at waveugide - set PML_opt as 1 for PML waveguide else PEC at
% boundary
PML_opt = 0; % If one wants PML inside the waveguide
PMLWG = 20;
S_WG = 20;
[Sigmax,isPEC] = PML_WG(Sigmax,WG0,n,a,PMLWG,S_WG,isPEC,PML_opt);

isPEC = logical(isPEC);
Sigmay(isPEC) = S_PEC; %PEC
Sigmax(isPEC) = S_PEC; %PEC

% Magnetic conductivity
Sigmay_star = Sigmay.*mur./er;
Sigmax_star = Sigmax.*mur./er;
Sigmax_star(isPEC) = 0;
Sigmay_star(isPEC) = 0;

% Init fields
Ezx = 0.*I_region;
Ezy  = Ezx;
Hx = Ezx;
Hy = Ezx;

% time and pulse
time = 1.5e3;
t0 = 1/(pi*0.625e9);
count = 1;

% Multiplication for fields
% % for H field
A = (1-Sigmay_star.*dt./(2.*mur))./(1+Sigmay_star.*dt./(2.*mur)); %Hx
B = (dt/dx)./(mur.*(1+Sigmay_star.*dt./(2.*mur))); %Hx
C = (1-Sigmax_star.*dt./(2.*mur))./(1+Sigmax_star.*dt./(2.*mur)); %Hy
D = (dt/dx)./(mur.*(1+Sigmax_star.*dt./(2.*mur))); % Hy

% for E field
E = (1-(Sigmax*dt)./(2.*er))./(1+(Sigmax*dt)./(2*er));
F = dt./(er.*dx)./(1+(Sigmax*dt)./(2*er));
G = (1-(Sigmay*dt)./(2.*er))./(1+(Sigmay*dt)./(2*er));
H = dt./(er.*dx)./(1+(Sigmay*dt)./(2*er));

% Deciding on where to place ficticious surface
Surf = ceil(1*wl/dx);
FF_surf = [PML+10,x_cell-PML-10;PML+10,y_cell-10-PML];
FF_surf = round(FF_surf);
x1 = FF_surf(1,1);
x2 = FF_surf(1,2);
y1 = FF_surf(2,1);
y2 = FF_surf(2,2);
Ezy2_sum = 0;
count = 1;
%Preallocating
sanity_check = zeros(1,time);
sctmp = zeros(1,time);
absSvec = zeros(1,time);
% E Field At the artificial Surface
Ezx1 = zeros(time,FF_surf(1,2) - FF_surf(1,1) + 1);
Ezx2 = zeros(time,FF_surf(1,2) - FF_surf(1,1) + 1);
Ezy1= zeros(time,FF_surf(2,2) - FF_surf(2,1) + 1); % Near the waveuide
Ezy2 = zeros(time,FF_surf(2,2) - FF_surf(2,1) + 1); % At the opening of antenna
% H Field At the artificial Surface
Hx1 = zeros(time,FF_surf(1,2) - FF_surf(1,1) + 1);
Hx2 = zeros(time,FF_surf(1,2) - FF_surf(1,1) + 1);
Hy1 = zeros(time,FF_surf(2,2) - FF_surf(2,1) + 1);
Hy2 = zeros(time,FF_surf(2,2) - FF_surf(2,1) + 1);

for t = 1:time
    % Solving Maxwells FDTD Equations - Yee Staggered
    Hy(1:x_cell-1,1:y_cell-1) = C(1:x_cell-1,1:y_cell-1).*Hy(1:x_cell-1,1:y_cell-1) ...
        + D(1:x_cell-1,1:y_cell-1).*(Ezx(2:x_cell,1:y_cell-1) + Ezy(2:x_cell,1:y_cell-1) ...
        -((Ezx(1:x_cell-1,1:y_cell-1) + Ezy(1:x_cell-1,1:y_cell-1))));
    
    Hx(1:x_cell-1,1:y_cell-1) = A(1:x_cell-1,1:y_cell-1).*Hx(1:x_cell-1,1:y_cell-1) ...
        +B(1:x_cell-1,1:y_cell-1).*(Ezx(1:x_cell-1,1:y_cell-1) + Ezy(1:x_cell-1,1:y_cell-1) ...
        -((Ezx(1:x_cell-1,2:y_cell) + Ezy(1:x_cell-1,2:y_cell))));
    Ezx(2:x_cell,2:y_cell) = E(2:x_cell,2:y_cell).*Ezx(2:x_cell,2:y_cell) + F(2:x_cell,2:y_cell) ...
        .*(Hy(2:x_cell,2:y_cell) - Hy(1:x_cell-1,2:y_cell));
    
    Ezy(2:x_cell,2:y_cell) = G(2:x_cell,2:y_cell).*Ezy(2:x_cell,2:y_cell) - H(2:x_cell,2:y_cell) ...
        .*(Hx(2:x_cell,2:y_cell) - Hx(2:x_cell,1:y_cell-1));
    
    % Hard source
    source= 1.484*exp(-((t*dt-3*t0)/t0)^2)*sin(2*pi*1e9*dt*t);
    sctmp(t) = source;
    Ezx(S0(1),S0(2)) = source;
    Ezy(S0(1),S0(2)) = source;
    
    Ez = Ezx+Ezy;
    Ezhist{count} = Ez;
    count = count+1;
    Ezx1(t,:) = Ez(x1:x2,y1); 
    Ezx2(t,:) = Ez(x1:x2,y2);
    Ezy1(t,:) = Ez(x1,y1:y2); % Near the waveuide
    Ezy2(t,:) = Ez(x2,y1:y2); % At the opening of antenna
    % H Field At the Surface
    Hx1(t,:) = (Hx(x1:x2,y1) + Hx(x1:x2,y1-1))/2;
    Hx2(t,:) = (Hx(x1:x2,y2) + Hx(x1:x2,y2-1))/2;
    Hy1(t,:) = (Hy(x1,y1:y2) + Hy(x1-1,y1:y2))/2;
    Hy2(t,:) = (Hy(x2,y1:y2) + Hy(x2,y1:y2))/2;
    
    % Sanity check -  Looking at E-field in one cell
    sanity_check(t) = Ez(340,y_cell/2);
    
    % Poynting Vector to check conservation
    Sex = -Ez.*Hy;
    Sey = Ez.*Hx;
    absScell = sqrt(Sex.^2 + Sey.^2);
    absS = sum(sum(absScell));
    absSvec(t) = absS;
end
plot(sanity_check)
hold on
plot(sctmp)
% fixing H on the surface to match E field time
[~,b] = size(Hx1);
Hx1_tmp = zeros(t,b);
Hx2_tmp = zeros(t,b);
for p = 1:b
    A = Hx1(:,p);
    A = A.';
    Hx1_tmp(2:end,p) = mean([A(1:end-1);A(2:end)]);
    A = Hx2(:,p);
    A = A.';
    Hx2_tmp(2:end,p) = mean([A(1:end-1);A(2:end)]);
end

[~,b] = size(Hy1);
Hy1_tmp = zeros(t,b);
Hy2_tmp = zeros(t,b);
for p = 1:b
    A = Hy1(:,p);
    A = A.';
    Hy1_tmp(2:end,p) = mean([A(1:end-1);A(2:end)]);
    A = Hy2(:,p);
    A = A.';
    Hy2_tmp(2:end,p) = mean([A(1:end-1);A(2:end)]);
end
%% Energy Conservation
plot(absSvec./max(abs(absSvec)))
hold on
plot([70,70],[-1,1])
hold on
plot([475,475],[-1,1])
legend('Poynting Vector absolute value','Source stop exciting','PML')
set(gca,'Fontsize',20)
%% Convergance by checking E-field
ni = [400,450,500,800,1e3];
plotting{1} = 'm';
plotting{2} = 'r';
plotting{3} = 'g';
plotting{4} = 'y';
plotting{5} = 'b--';
for k = 1:length(ni)
    clear Ezz
    clear Ez_fft
    clear Ez_abs
    n = ni(k);
    Ezz = Ezy2(1:n,:);
    [Ez_fft,s_fft] = get_fft(Ezz,1.6e9);
    Ez_abs = abs(Ez_fft);
    plot(1:s_fft,Ez_abs,plotting{k},'linewidth',1.2)
    hold on
end
ylabel('Relative Electric Field Strength')
xlabel('Cell number')
legend('t = 400','t = 450','t = 500','t = 800','t = 1000')
set(gca,'FontSize',20)
%% Compare near field cst and nearfield matlab FDTD
fid= fopen('nearfield_16ghz_cst.txt');
[a,b] = cst_txt(fid);
clear Ezz
clear Ez_fft
clear Ez_abs
[Ez_fft,s_fft] = get_fft(Ezy2,1.6e9);
Ez_abs = abs(Ez_fft);
plot(1e3*dx*(-s_fft/2+1:s_fft/2),Ez_abs/(max(Ez_abs)),'r','linewidth',1.5)
hold on
plot(a,b/max(b),'b--','linewidth',1.5)
ylabel('Normalized E-field Strength,[V/V]')
xlabel('Distance [mm]')
legend('Matlab Simulation','CST Simulation')
set(gca,'fontsize',20)
% print -depsc comparison.eps
%% The phase at the artificial surface - Comparison CST and FDTD
fid= fopen('nearfield_phase_16.txt');
[a,b] = cst_txt(fid);
b = deg2rad(b);
b = unwrap(b);
phase = unwrap((angle(Ez_fft)));
diff = max(phase)-max(b);
plot(1e3*dx*(-s_fft/2+1:s_fft/2),phase-diff,'r','linewidth',1.5)
hold on
plot(a,b,'b--','linewidth',1.5)
legend('Matlab Simulation','CST Simulation')
xlabel('Distance [mm]')
ylabel('Phase [rad]')
set(gca,'FontSize',20)
% print -depsc comparison_phase.eps

%% Compare E-field between no lens and with lens
load phase_nl.mat
load phase_lens.mat

plot(1e3*dx*(-s_fft/2+1:s_fft/2),Ez_nl,'b','linewidth',1.5)
hold on
plot(1e3*dx*(-s_fft/2+1:s_fft/2),Ez_lens,'r','linewidth',1.5)

legend('Field without lens','Field with lens')
xlabel('Distance [mm]')
ylabel('Relative Electric Field Strength')
set(gca,'FontSize',20)
% print -depsc E_comp_lens.eps
%% Compare phase between no lens and with lens
diff = max(phase_nl)-max(phase_lens);

plot(1e3*dx*(-s_fft/2+1:s_fft/2),phase_nl - diff,'b','linewidth',1.5)
hold on
plot(1e3*dx*(-s_fft/2+1:s_fft/2),phase_lens,'r','linewidth',1.5)

legend('Phase without lens','Phase with lens')
xlabel('Distance [mm]')
ylabel('Phase [rad]')
set(gca,'FontSize',20)
% print -depsc phase_comp_lens.eps

%% NF -> FF Transformation
iso_dbi = -12.11-2.15; % Gaussian pulse gives directivity of -13.6, we use this as a refrence
fid= fopen('ff_1_6GHz.txt');
[aa,bb] = cst_txt_ff(fid);
% aa = linspace(0,2*pi,length(aa));
D_horn = l_rwg + R;
f0 = 1.6e9;
wl0 = c0/f0;
k = 2*pi/wl0; % Wave number
r_FF = 4*(D_horn*dx)^2/(wl0*f/f0); % Farfield distance in meters
L1 = [x1,x2; y1,y1];
L2 = [x2,x2; y1,y2];
L3 = [x1,x2; y2,y2];
L4 = [x1,x1;y1,y2];
% rectangle('Position',[x1,y1,x2-x1,y2-y1])
Hx1 = Hx1_tmp;
Hx2 = Hx2_tmp;
Hy1 = Hy1_tmp;
Hy2 = Hy2_tmp;
[Ez_L1,~] = get_fft(Ezx1,f0);
[Ez_L3,~] = get_fft(Ezx2,f0);
[Ez_L4,~] = get_fft(Ezy1,f0);
[Ez_L2,~] = get_fft(Ezy2,f0);

[Hx_L1,~] = get_fft(Hx1,f0);
[Hx_L3,~] = get_fft(Hx2,f0);
[Hy_L4,~] = get_fft(Hy1,f0);
[Hy_L2,~] = get_fft(Hy2,f0);

phi_L1 = 0;
phi_L4 = pi/2;
phi_L2 = atan((0:(y2-y1))./(x2-x1)).';
phi_L3 = atan((y2-y1)./((x2-x1):-1:0)).';

L1_dist = 0:(x2-x1);
L1_dist = L1_dist.'*dx;

L2_dist = ((x2-x1)^2 + (0:(y2-y1)).^2).^0.5;
L2_dist = L2_dist.'*dx;

L3_dist = ((y2-y1)^2 + ((x2-x1):-1:0).^2).^0.5;
L3_dist = L3_dist.'*dx;

L4_dist = (y2-y1):-1:0;
L4_dist = L4_dist.'*dx;

K1 = 2*pi*f0*mur0;
K2 = -exp(-1j*k*r_FF)/sqrt(r_FF)*sqrt(1j)/(sqrt(8*pi*k));
ang = linspace(0,2*pi,100);
for q = 1:length(ang)
    PHI = ang(q);
    psi_L1 = phi_L1-PHI;
    Eff_L1 = (K1.*Hx_L1 + k.*Ez_L1.*sin(PHI)).*exp(1j*k.*L1_dist.*cos(psi_L1))*dx;
    Ez_L1_int(q) = sum(Eff_L1);
    
    psi_L2 = phi_L2-PHI;
    Eff_L2 = (K1.*Hy_L2 - k.*Ez_L2.*cos(PHI)).*exp(1j*k.*L2_dist.*cos(psi_L2))*dx;
    Ez_L2_int(q) = sum(Eff_L2);
    
    psi_L3 = phi_L3-PHI;
    Eff_L3 = (-K1.*Hx_L3 - k.*Ez_L3.*sin(PHI)).*exp(1j*k.*L3_dist.*cos(psi_L3))*dx;
    Ez_L3_int(q) = sum(Eff_L3);
    
    psi_L4 = phi_L4-PHI;
    Eff_L4 = (-K1.*Hy_L4 + k.*Ez_L4.*cos(PHI)).*exp(1j*k.*L4_dist.*cos(psi_L4))*dx;
    Ez_L4_int(q) = sum(Eff_L4);     
end
INT = Ez_L1_int + Ez_L2_int + Ez_L3_int + Ez_L4_int;
INT = K2.*INT;
INT_abs = abs(INT);
FF_db = 10*log10(INT_abs);

plot(rad2deg(ang),FF_db-iso_dbi,'r','linewidth',1.5)
hold on
plot(aa,bb,'b','linewidth',1.5)

xlabel('\phi [Degrees]')
ylabel('Directivity [dBi]')
legend('Farfield from Matlab','Farfield from CST')
set(gca,'Fontsize',20)
% print -depsc ff_comparison.eps
%% Compare Farfield of lens and no lens
iso_dbi = -12.11-2.15; % Gaussian pulse gives directivity of -13.6, we use this as a refrence
load FF_lens.mat
load FF_nl.mat
ang = linspace(0,2*pi,100);
plot(rad2deg(ang),FF_nl -iso_dbi,'b','linewidth',1.5)
hold on
plot(rad2deg(ang),FF_lens -iso_dbi,'r','linewidth',1.5)
xlabel('\phi [Degrees]')
ylabel('Directivity [dBi]')
legend('Farfield without lens','Farfield with lens')
set(gca,'FontSize',20)
% print -depsc FF_lens.eps
%% Compare farfield but for wider aperture
iso_dbi = -12.11-2.15; % Gaussian pulse gives directivity of -13.6, we use this as a refrence
load FF_lens_wider.mat
load FF_nl_wider.mat
ang = linspace(0,2*pi,100);
plot(rad2deg(ang),FF_nl -iso_dbi,'b','linewidth',1.5)
hold on
plot(rad2deg(ang),FF_lens -iso_dbi,'r','linewidth',1.5)
xlabel('\phi [Degrees]')
ylabel('Directivity [dBi]')
legend('Farfield without lens','Farfield with lens')
set(gca,'FontSize',20)
% print -depsc FF_lens_wider.eps
%% plotting the Ez - imagesc timeplot of the Ez field
for ss = 1:t
    if mod(ss,15) == 0 || ss == 1
        imagesc(Ezhist{ss})
        caxis([-0.005,0.005])
        colorbar
        pause(0.0005)
    end
end
%% Comparing PML and PEC - Comparison of the E field with a perfect mathcing layer and perfect metal cond 
% The .mat files aren't included
count2 = 50;
load pml_ez.mat
load PEC_ez.mat
for ss = 1:count2
    subplot(1,2,1)
    imagesc(Ez_pml{ss})
    caxis([-0.001,0.001])
    title('PML BC')
    xlabel('Cell y-axis')
    ylabel('Cell x-axis')
    c = colorbar;
    title(c, '\fontsize{15} E_z')
    set(gca,'fontsize',20)
    subplot(1,2,2)
    imagesc(Ez_PEC{ss})
    caxis([-0.001,0.001])
    title('PEC BC')
    xlabel('Cell y-axis')
    ylabel('Cell x-axis')
    set(gca,'fontsize',20)
    c = colorbar;
    title(c,'\fontsize{15} E_z')
    pause(0.05)
end
