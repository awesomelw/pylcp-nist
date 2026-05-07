%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   This code plots kvec and pol of a rotated gaussian beam, in cartesean
%   coordiates
%
%   Edited by Chad 07/20/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

% Parameters to change
wb=200;         %beam waist (nm);
beta = 1;       % Beam intensity in focus
kvec=[1,1,1]/sqrt(3)*(2*pi/461); % k-vector, normalized by 2pi/lambda(nm)
pol=[1,+1i,0]/sqrt(2);   % keeping track of LHP (assuming z-axis reference)
r0=[100,100,0];     % Beam offset

% Define grid for evaluationg Gaussian function
xa=[-600:200:600];
ya=[-600:200:600];
za=[-600:200:600];
[x,y,z]=meshgrid(xa,ya,za);
r=[x(:),y(:),z(:)];

% Calculate for both x and y polarizations
[I,k,P]=GaussCalc(kvec,pol,r0,beta,wb,r);

%% Check results
Ip=I; % scale factor for plotting
Imax=0.5*max(Ip);
Ip(Ip>Imax)=Imax; % Threshold intensity for easier visualization

figure
% scatter3(r(:,1),r(:,2),r(:,3),I*100,'filled','k');hold on
quiver3(r(:,1),r(:,2),r(:,3),Ip.*k(:,1),Ip.*k(:,2),Ip.*k(:,3),'k');hold on; % k-vector
quiver3(r(:,1),r(:,2),r(:,3),Ip.*real(P(:,1)),Ip.*real(P(:,2)),Ip.*real(P(:,3)),'r');   % x-pol
quiver3(r(:,1),r(:,2),r(:,3),Ip.*imag(P(:,1)),Ip.*imag(P(:,2)),Ip.*imag(P(:,3)),'b');   % y-pol
axis equal;grid on;xlabel('x');ylabel('y');zlabel('z')

% Check to make sure all (normalized) vectors are orthogonal
for n=1:size(k,1)
    Dot(n)=abs(dot(k(n,:)/norm(k(n,:)),real(P(n,:))/norm(real(P(n,:)))))+...
        abs(dot(k(n,:)/norm(k(n,:)),imag(P(n,:))/norm(imag(P(n,:)))))+...
        abs(dot(real(P(n,:))/norm(real(P(n,:))),imag(P(n,:))/norm(imag(P(n,:)))));
end
maxdotproduct=max(Dot)






%% Function to call
function [I,k,P]=GaussCalc(kvec,pol,r0,beta,wb,r)

% k-vect
kmag=norm(kvec);
khat=kvec/kmag;
wavelength=2*pi/kmag;

% Translation
Rt=r-r0;

% Global Rotation matrix
th = acos(khat(3));
ph = atan2(khat(2),khat(1));
rz=[cos(ph),-sin(ph),0;sin(ph),cos(ph),0;0,0,1];
ry=[cos(th),0,sin(th);0,1,0;-sin(th),0,cos(th)];
% rmat=rz*ry;           %%% OLD definition of rotation matrix (works for k, but doesn't work for polarization)
rmat=rz*ry*inv(rz);     %%% NEW definition of rotation matrix (works for k, preserves polarizations well)
rmatinv=inv(rmat);
for n=1:size(Rt,1)
    Rtr(n,:)=rmatinv*transpose(Rt(n,:));    % I need to do an inverse rotatoin on the positions
end

% Save translated and rotated positions
x=Rtr(:,1);
y=Rtr(:,2);
z=Rtr(:,3);
rho=sqrt(x.^2+y.^2);

% Gaussian sub functions
zr=pi*wb^2/wavelength;
w=wb*sqrt(1+(z/zr).^2);

%% Intensity
I=beta*(wb./w).^2.*exp(-2*rho.^2./w.^2);

%% k-vector (from L. Allen, M.J. Padgett / Optics Communications 184 (2000) 67-71)
kr=rho.*z./(zr^2+z.^2);   % kr = kr/kz from reference ^
kt=0;  % Assuming w>lambda, so kt/kz ~= 0
kz=ones(size(rho));       % Normalizing kz/kz = 1;

% changing to cartesian coordinates
kx=(kr.*x+kt.*y)./sqrt(x.^2+y.^2+1e-100); %Added a small number to prevent 0/0 = NaN condition
ky=(kr.*y+kt.*x)./sqrt(x.^2+y.^2+1e-100);

% Normalize k-vector for plotting
kxn=kx./sqrt(kx.^2+ky.^2+kz.^2);
kyn=ky./sqrt(kx.^2+ky.^2+kz.^2);
kzn=kz./sqrt(kx.^2+ky.^2+kz.^2);

% Rotate k-vector and polarizations
for n=1:size(Rt,1)   
    
    % Local rotation based on local k-vector
    thn = acos(kzn(n));
    phn = atan2(kyn(n),kxn(n));
    rzn=[cos(phn),-sin(phn),0;sin(phn),cos(phn),0;0,0,1];
    ryn=[cos(thn),0,sin(thn);0,1,0;-sin(thn),0,cos(thn)];
%     rmatn=rz*ry;              %%% Old rotation matrix (doesn't work well)
    rmatn=rzn*ryn*inv(rzn);     %%% New rotation matrix for local polarization
    
    % Local rotate then global rotate polarization vectors
    P(n,:)=rmat*rmatn*transpose(pol);
    
    % Global rotate k-vect
    Kn(n,:)=rmat*[kxn(n);kyn(n);kzn(n)];
end

k=Kn*kmag;  %Renormalize rotated k-vector with kmagnitude

end