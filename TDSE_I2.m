%% TDSE solver using split step method for wavepacket dynamics in I2
% --------------------------------------------------------------------
%  
%   Adi Natan (natan@stanford.edu)
%
% --------------------------------------------------------------------

% Turn on (1) and off (0) plotting of the wavepacket propagation
plot_flag=0;

%constants & conversions
h_bar=1; i1=sqrt(-1);
AU_TO_EV=27.2113961;
NM_TO_AU=45.56335;
AU_TO_CM=2.194746e5;
AU_TO_ANG=0.529177249;
AU_TO_FS=0.0241889;
AMU_TO_AU=1.8228885e3;
WCM2_TO_AU=3.5094452*1e16;
conv_r= AU_TO_ANG;
conv_en=AU_TO_CM;
conv_t= AU_TO_FS;

% set up R grid (in Bohr)
N=2^11; % grid points 
Rmin = eps;  Rmax = 40; % box size in au
R = linspace(Rmin,Rmax,N); % in au
dR = R(2)-R(1);

r=R*conv_r; % R grid in Ang

% mass of I2
mass = (126.904473/2)*AMU_TO_AU;

% set up k grid
dk = 2*pi/(N*dR);
k =[(0:N/2-1),(-N/2:-1)]*dk;

% set up kinetic energy and k-space propagator
kin= (0.5*(k.^2)/mass).';

% set up the potential curves od I2:
% X (ground) state of I2
De_X= 12440/AU_TO_CM;
be_X= 1.875*AU_TO_ANG;
Re_X= 2.656/AU_TO_ANG;
V1= De_X*(1-exp(-be_X*(R-Re_X))).^2;
V1=V1(:);

% B (bound excited) state of I2
De_B=  5169/AU_TO_CM;
be_B= 1.696*AU_TO_ANG;
Re_B= 3.025/AU_TO_ANG;
Te_B= 15771/AU_TO_CM;
V2= De_B*(1-exp(-be_B*(R-Re_B))).^2 + Te_B;
V2=V2(:);

% C (dissociative excited) state of I2
De_C=.0543;
alpha_C= -3.6060;
be_C=  1.6041;
Re_C=3.5161;
Te_C= 522.5292/AU_TO_CM;
V3= De_C*(1-alpha_C*exp(-be_C*(R-Re_C))).^2 + Te_C;
V3=V3(:);

%% get initial psi, use eigenstate of ground state via Fourier Grid Method
kk=pi/dR;
% kinetic energy
[nx, ny]=meshgrid(1:N);
T=-2*kk^2*(-1).^(nx-ny)./(pi*(nx-ny)).^2;
T(1:N+1:end)=0; % zeros inf on diagonal
T=-(0.5*h_bar^2/mass)*(T-eye(N)*kk^2/3);

% get the ground state (GS) of I2
[evec1,evl1]=eig(T+diag(V1)); % eigen-states of X-state
norm1= 1/sqrt(dot(evec1(:,1), evec1(:,1)));
psi_GS= norm1*evec1(:,1);    
E_GS=evl1(1,1);
%% time propagation parameters 
% laser pulse parameters
I0     =  1e11;            % peak power in W/cm^2
FWHM   =  30.0/AU_TO_FS;   % pulse FWHM duration
lambda =  520;             % pulse wavelength in nm
w0     =  NM_TO_AU/lambda; % pulse frequency
t0     =  3*FWHM;          % pulse envelope center

% Pulse intensities in atomic units
E0 = sqrt(I0/WCM2_TO_AU);    % the peak field amplitude in a.u. (1st pulse)

% Transition dipole in atomic units
mu12 = 0.181;
mu13 = 0.181;   

% time grid parameters
dt       = 0.1/AU_TO_FS; 
t_sample = 15/AU_TO_FS;
t_end    = 1000/AU_TO_FS;
n_steps  = round(t_sample/dt);
n_sample = round(t_end/t_sample);
t_end    = n_sample*t_sample;
Nplt     = 300;
t        = 0:dt:t_end; % time grid

% The laser pulse
E=E0*(exp(-2*log(2)/FWHM^2*(t-t0).^2).*exp(-1i*w0*(t-t0)));

%%  propagate !
psi1=    psi_GS;          % The GS unperturbed wavefunction
psi2=    zeros(N, 1);
psi3=    zeros(N, 1);

if plot_flag
    hplt=plot(r,abs(psi1).^2,...
        r,V2(end)+abs(psi2).^2, ...
        r,V3(end)+abs(psi3).^2, ...
        r,V1,r,V2,r,V3);
    ylabel('Energy (au)');
    xlabel('r [A]')
    ylim([0 0.15])
    xlim([2 8])

    drawnow;
end

Uk2=exp(-i1*kin*abs(dt)/2);
UV12=exp(-i1*(V1+V2)*abs(dt)/2);
UV13=exp(-i1*(V1+V3)*abs(dt)/2);
UV23=exp(-i1*(V2+V3)*abs(dt)/2);

n=0;

textprogressbar('Generate charge density rho(r,tau), solving TDSE for diatomic I2:'); tic

for idx=1:numel(t)

    [psi1, psi2]=prop_dt2(dt,mu12*E(idx),psi1,psi2,V1,V2,Uk2,UV12);
    [psi1, psi3]=prop_dt2(dt,mu13*E(idx),psi1,psi3,V1,V3,Uk2,UV13);

    nr=sqrt(sum(abs(psi1).^2 + abs(psi2).^2+ abs(psi3).^2));
    psi1=psi1./nr;
    psi2=psi2./nr;
    psi3=psi3./nr;

    if mod(idx-1,Nplt)==0
        if plot_flag
            % just for visulaization, we scale the ground and excited
            % states and dress them on their potential curves,
            set(hplt(1),'YData', 0.25*abs(psi1))
            set(hplt(2),'YData', V2(end)+5*abs(psi2).^2);
            set(hplt(3),'YData', V3(end)+5*abs(psi3).^2);
            title(['t= ' num2str(t(idx)*conv_t) ' fs']);
            legend('X - Ground State (x0.25)','B - Excited State (x5)','C - Excited State (x5)')
            drawnow;
        end
        n=n+1;
        rho1(:,n)=abs(psi1).^2;
        rho2(:,n)=abs(psi2).^2;
        rho3(:,n)=abs(psi3).^2;
        time(n)=t(idx)*conv_t;

    end
    textprogressbar(idx./numel(t)*100);

end
textprogressbar(['done! ' '(' num2str(toc) ' sec)']);

rho=rho1+rho2+2*rho3; % adding more to the diss. wavepacet for visibility

%% plot
figure(100)
imagesc(time,r,rho-rho(:,1));colorbar
caxis([-3e-3 3e-3]);
 
set(gca,'Ydir','normal')
xlabel('\tau [fs]');
ylabel('R [A]');
title('\Delta\rho(R,\tau)')
ylim([2 20])

%% save
t=time-time(4); % set t=0 for pulse center;
save('TDSE_I2.mat','rho', 't','r');

function [psi1, psi2]=prop_dt2(dt,Om,psi1,psi2,V1,V2,Uk2,UV12)
% propagate two wavepackets via split step fourier, see D. Tannor's book
% Introduction to Quantum Mechanics: A Time-Dependent Perspective, 11.7.1
    i1=sqrt(-1);
    D=sqrt((V2-V1).^2+abs(4*abs(Om).^2));
    cosD=cos(D*dt/2);
    sinDovD=sin(D*dt/2)./D;
    
    psi1=ifft(Uk2.*fft(psi1));
    psi2=ifft(Uk2.*fft(psi2));
    
    psi1tmp=UV12.*((cosD+i1*(V2-V1).*sinDovD).*psi1-2*i1*Om.*sinDovD.*psi2);
    psi2tmp=UV12.*((cosD+i1*(V1-V2).*sinDovD).*psi2-2*i1*Om.*sinDovD.*psi1);
    
    psi1=ifft(Uk2.*fft(psi1tmp));
    psi2=ifft(Uk2.*fft(psi2tmp));
end   
 
