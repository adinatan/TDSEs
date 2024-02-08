%% TDSE solver using split step method for wavepacket pump-dump perturbative dynamics in I2
% --------------------------------------------------------------------
%
%   Adi Natan (natan@stanford.edu)
%   Feb 7, 2024
% --------------------------------------------------------------------
ccc
% Turn on (1) and off (0) plotting of the wavepacket propagation
plot_flag=1;

%constants & conversions
h_bar=1;
i1=sqrt(-1);
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
N=2^9; % grid points
Rmin = 2;  Rmax = 11; % box size in au
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

%% set up the potential curves od I2:
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

%% get initial psi, use eigenstate of ground state via Fourier Grid Method
kk=pi/dR;
% kinetic energy
[nx, ny]=meshgrid(1:N);
T=-2*kk^2*(-1).^(nx-ny)./(pi*(nx-ny)).^2;
T(1:N+1:end)=0; % zeros inf on diagonal
T=-(0.5*h_bar^2/mass)*(T-eye(N)*kk^2/3);

% get the ground state (GS) of I2
[evec1,evl1]=eig(T+diag(V1)); % eigen-states of X-state
Evl1=diag(evl1); % eigen energies of V1

% since I2 was at ~400K we actually need to coisder the following vibrarions dist frac.
% v=0 53.58%
% v=1 24.87%
% v=2 11.55%
% v=3 5.36%
% v=4 2.49%
% so we should consider later to loop over each and add
% incoherently the result with the assigned weights in vf:
%vf=[0.5358 0.2487 0.1155 0.0536 0.0249];
 
v_init=0; % this is to test for the cold GS v=0,

psi_GS= (evec1(:,v_init+1))';
psi_GS= psi_GS/sqrt(dot(psi_GS,psi_GS))';
E_GS=Evl1(v_init+1);

%% time propagation parameters
% laser pulse parameters
I0_1     =  1e11;            % peak power in W/cm^2
I0_2     =  1e11;            % peak power in W/cm^2

FWHM   =  30.0/AU_TO_FS;   % pulse FWHM duration
lambda_1 =  520;             % pulse wavelength in nm
lambda_2 =  800;             % pulse wavelength in nm

w0_1     =  NM_TO_AU/lambda_1; % pulse frequency
w0_2     =  NM_TO_AU/lambda_2; % pulse frequency

t0_1     =  2*FWHM;   % pulse 1 envelope center
te =  450/AU_TO_FS;   % Time delay between pump and dump pulses
t0_2=  t0_1 + te;     % pulse 2 envelope center

phi_2= 0; % phase of w0_2, possibly run over several phases and avg.

% Pulse intensities in atomic units
E0_1 = sqrt(I0_1/WCM2_TO_AU);    % the peak field amplitude in a.u. (1st pulse)
E0_2 = sqrt(I0_2/WCM2_TO_AU);    % the peak field amplitude in a.u. (2nd pulse)


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
Nplt     = 100;
t        = 0:dt:t_end; % time grid


% The laser pulses
pulse1= E0_1*exp(-2*log(2)/FWHM^2*(t-t0_1).^2).*exp(-1i*w0_1*(t-t0_1));
pulse2= E0_2*exp(-2*log(2)/FWHM^2*(t-t0_2).^2).*exp( 1i*w0_2*(t-t0_2))*exp(1i*phi_2);


%% set up propagators
%u_kin= exp(-1i*dt*kin);
u_kin_half= exp(-1i*0.5*dt*kin);
u_pot1= exp(-1i*dt*(V1));
u_pot2= exp(-1i*dt*(V2));
%u_pot3= exp(-1i*dt*(V3));

psi1=    psi_GS.';          % unperturbed wavefunction
psi2=    zeros(N, 1);     % First order wavefunction
psi21=   zeros(N, 1);     % second order wavefunction
%psi3=  zeros(N, 1);     % third order wavefunction

% only plot excited wavepackets
if plot_flag
    hplt=plot(r,V2(end)+abs(psi2).^2, ...
        r,E_GS+w0_1-w0_2+abs(psi21).^2,...
        r,V1,r,V2);

    ylabel('Energy (au)');
    xlabel('r [A]')
    ylim([0 0.15])
    xlim([2.1 5])

    drawnow;
end

n=0;
textprogressbar('Generate charge density rho(r,tau), solving TDSE for diatomic I2:'); tic

rho1=zeros(N,ceil(numel(t)/Nplt));
rho2=rho1;
rho21=rho2;

for idx=1:numel(t)

    %while (tot_steps<max_steps)
    % propagate wavefunctions 1/2 step forward
    U= u_kin_half;
    psi2=    ifft(U.*fft(psi2));
    psi21=   ifft(U.*fft(psi21));

    % potential step
    psi1=    psi1*exp(-1i*dt*E_GS);
    psi2=    u_pot2.*psi2;
    psi21=   u_pot1.*psi21;


    psi2=    ifft(U.*fft(psi2));
    psi21=   ifft(U.*fft(psi21));

    % Add contribution due to field interactions
    psi2=    psi2    + 1i*dt*(pulse1(idx)*mu12.*psi1);
    psi21=   psi21   + 1i*dt*(pulse2(idx)*mu12.*psi2);

    % not relevant because we use a perturbative approach
    % nr=sqrt(sum(abs(psi1).^2 + abs(psi2).^2+ abs(psi21).^2));
    % psi1=psi1./nr;
    % psi2=psi2./nr;
    % psi21=psi21./nr;

    % track populations
    pop2= dot(psi2,psi2);
    pop_2(idx)= pop2;
    r_1(idx)= real(dot(psi2, (r.').*psi2)/pop2);

    pop21= dot(psi21,psi21);
    pop_21(idx)= pop21;
    r_21(idx)= real(dot(psi21, (r.').*psi21)/pop21);


    if mod(idx-1,Nplt)==0
        if plot_flag

            % just for visulaization, we scale the ground and excited
            % states and dress them on their potential curves,
            set(hplt(1),'YData', V2(end)+10*abs(psi2).^2);
            set(hplt(2),'YData', E_GS+w0_1-w0_2+5e2*abs(psi21).^2);
            legend('B Excited state (x10)','Raman dump (x500)','Location','southeast')
            drawnow;

            title(  ['(t-t_1)=' num2str(conv_t*(t(idx)-t0_1)) ...
                ' , (t-t_2)= ' num2str(conv_t*(t(idx)-t0_2))]);
            drawnow

        end
        n=n+1;
        rho1(:,n)=abs(psi1).^2;
        rho2(:,n)=abs(psi2).^2;
        rho21(:,n)=abs(psi21).^2;
        time(n)=t(idx)*conv_t;
    end
    %textprogressbar(idx./numel(t)*100);

end
textprogressbar(['done! ' '(' num2str(toc) ' sec)']);



%% Plot pulses and populations
figure(2);
tiledlayout("flow")
nexttile
plot(t*conv_t, abs(pulse1+pulse2).^2);
ylabel('Field intensity');
xlabel('\tau [fs]');

nexttile
plot(t*conv_t, pop_2);
ylabel('B state population');
xlabel('\tau [fs]');

nexttile
plot(t*conv_t, pop_21);
ylabel('Raman population');
xlabel('\tau [fs]');


nexttile
imagesc(time,r,rho2);colorbar
%caxis([-3e-3 3e-3]);
colorbar
set(gca,'Ydir','normal')
xlabel('\tau [fs]');
ylabel('R [A]');
title('\rho_2(R,\tau)')
ylim([2 5])


nexttile
imagesc(time,r,rho21);colorbar
%caxis([-3e-3 3e-3]);
colorbar
set(gca,'Ydir','normal')
xlabel('\tau [fs]');
ylabel('R [A]');
title('\rho_{21}(R,\tau)')
ylim([2 5])