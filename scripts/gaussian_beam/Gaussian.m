close all
clear
clc

% Parameters
lambda=461; %wavelength (nm);
w0=500;    %beam waist (nm);
Beta = 10;  % Beam intensity in focus

% Define grid for evaluationg Gaussian function
xa=[-1000:200:1000];
ya=[-1000:200:1000];
za=[-3000:1000:3000];
[x,y,z]=meshgrid(xa,ya,za);
r=sqrt(x.^2+y.^2);

% Gaussian sub functions
zr=pi*w0^2/lambda;
w=w0*sqrt(1+(z/zr).^2);

%% Intensity
I=Beta*(w0./w).^2.*exp(-2*r.^2./w.^2);

%% k-vector (from L. Allen, M.J. Padgett / Optics Communications 184 (2000) 67-71)
kr=r.*z./(zr^2+z.^2);   % kr = kr/kz from reference ^
kt=zeros(size(r));      % Assuming w>lambda, so kt/kz ~= 0
kz=ones(size(r));       % Normalizing kz/kz = 1;

% changing to cartesian coordinates
kx=(kr.*x+kt.*y)./sqrt(x.^2+y.^2+1e-100); %Added a small number to prevent 0/0 = NaN condition
ky=(kr.*y+kt.*x)./sqrt(x.^2+y.^2+1e-100);

% Normalize k-vector for plotting
kxn=kx./sqrt(kx.^2+ky.^2+kz.^2);
kyn=ky./sqrt(kx.^2+ky.^2+kz.^2);
kzn=kz./sqrt(kx.^2+ky.^2+kz.^2);


%% Plot
quiver3(x(:),y(:),z(:),I(:).*kxn(:),I(:).*kyn(:),I(:).*kzn(:),'k');
axis equal;view(-30,30);grid off
