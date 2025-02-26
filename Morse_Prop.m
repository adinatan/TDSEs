% Reference: Szczepan Chelkowski et al doi.org/10.1103/PhysRevLett.65.2355
clear all
close all
plot_flag = 'true';

% ================= conversions and constants ====================
% units conversions
au2amu = 1/1822.888; % m_e to amu
au2ang = 0.52917721092; % a_0 to angstrom
au2eV = 27.211385; % E_h to eV
au2as = 24.18884326505; % hbar/E_h to attosecond
au2Vpang = 51.4220652; % E_h/(e *a_0) to V/angstrom
au2D = 2.541746; % e *a_0 to Debye
J2eV = 1/1.6021766208e-19;
s2as = 1e18;
cm2ang = 1e8;

% constants in au
m_e = 1; % electron mass
e = 1;
hbar = 1;
k_e = 1;

a_0 = 1; % hbar^2/(k_e *m_e *e^2), Bohr radius
c = 137.035999139; % a_0*E_h/hbar = alpha *c
% constants in SI
k_e_SI = 8.9875517873681e9; % kg*m^3/(s^2 C^2)
c_SI = 2.99792458e8; % m/s

% ================= Laser Pulse ==========================
I = 1e13; % W/cm^2
E_M = sqrt(8*pi*k_e_SI/c_SI *I *1e4) *1e-10; % V/angstrom
E_M = E_M/au2Vpang; % au = E_h/(e*a_0)

m_H = 1.00794; % amu
m_F = 18.998403; % amu
m = m_H*m_F/(m_H +m_F); % amu
m = m/au2amu; % au, reduced mass

a = 1.1741/a_0; % au = 1/a_0
D = 6.125/au2eV; % au = E_h
B = hbar *a/sqrt(2*m*D); % dimless
omega_0 = 2*B*D/hbar; % au = E_h/hbar
omega_01 = omega_0 *(1-B); % au = E_h/hbar
cycle_01 = 2*pi/omega_01; % au = hbar/E_h

N = 8; % dimless
alpha_0 = 2.5; % dimless
alpha_F = 6.25; % dimless
S = 1.5*pi; S_0 = 1.5*pi; % dimless
d_1 = 0.786/(a_0*au2D); % au = e
p_01 = 0.097/au2D; % au = e*a_0

Q = S*hbar/(p_01 *E_M); % au = hbar/E_h
A = (2*(pi/4 - atan2(1,exp(alpha_0)))/alpha_0-1) /(1-sech(alpha_0));

t_0 = Q *(S/S_0 - 2*(sqrt(2) -sqrt(3/2)))/(A+1); % au = hbar/E_h
t_1 = Q *S/S_0 - A*t_0; % au = hbar/E_h
t_N = 2*Q*(sqrt(N+1)-sqrt(2))+t_1;
t_c = t_N + 0.02*t_0;

% time steps
dt = 0.001 *cycle_01; % au = hbar/E_h
t_max = 120 *cycle_01; % au = hbar/E_h
t = 0:dt:t_max;
Nt = numel(t);
%t_max = max(t); % au = hbar/E_h

omega = -B*((t-t_1).^2/(4*Q^2) +sqrt(2)*(t-t_1)/Q +3/2) +1; % dimless
omega = omega * omega_0; % au = E_h/hbar
omega(t < t_0) = omega_01; % 1/s
omega(t >= t_N) = omega(find(t < t_N,1,'last'));

U = ones(1,Nt);
U(t < t_0) = (1-sech(alpha_0))^(-1)*(sech(alpha_0*(t(t<t_0)-t_0)/t_0) -sech(alpha_0)); % dimless
U(t > t_c) = (1-sech(alpha_F))^(-1)*(sech(alpha_F*(t(t>t_c)-t_c)/t_0) -sech(alpha_F)); % dimless

%% ======================= 2 initial state ============================
% ================= Set Morse Potential =================
r_0 = 1.7329*a_0; % au = a_0
Nr = 2^10; % originally 2^11...
r_max = 65; % au = a_0, originally 130...
r = linspace(0.1,r_max,Nr); % au = a_0
dr = r(2) -r(1);
x = r - r_0; % au = a_0
dx = dr;

% Morse potential at time t = 0
V = D*(1 - exp(-a*x)).^2; % E_h

% ======================= ground state ===================
lam = 1/B; % dimless
z = 2*lam *exp(-a*x);

n = 0;
psi_0 = z.^(lam -n -1/2) .*exp(-z/2) .*laguerreL(n,2*lam-2*n -1,z);
psi_0 = psi_0/(sqrt(dr)*norm(psi_0));

%% ================= 3 propagation (Split Step) ============
dti_cycle_01 = round(cycle_01/(2*dt)); % index

dp = 2*pi*hbar/(Nr*dr); % au = hbar/a_0
p = ((0:Nr-1) -Nr/2) *dp; % au = hbar/a_0
T = p.^2/(2*m); % au = E_h
T = ifftshift(T);

dt_hbar = dt/hbar;
psi_1 = psi_0;

if plot_flag; figure; end

count = 0;
for ti = 1:Nt
    V = D*(1 - exp(-a*x)).^2 - x *d_1*E_M*U(ti)*cos(omega(ti)*t(ti));
    psi_1 = exp(-1i *dt_hbar/2 *V) .*psi_1;
    psi_1 = ifft(exp(-1i *dt_hbar *T) .*fft(psi_1));
    psi_1 = exp(-1i *dt_hbar/2 *V) .*psi_1;
    psi_1 = psi_1/(sqrt(dr)*norm(psi_1)); % normalize

    if rem(ti,dti_cycle_01/5) == 0
        count = count +1;
        psi(count,:) = psi_1;
        if plot_flag
            plot(r, abs(psi_1).^2, r, V,'LineWidth',2);
            text(6,2.5,['t = ' num2str(ti*dt/cycle_01, '%.1f') ,' cycles']);
            xlabel('R (a_0)');
            ylim([0,5]); xlim([0,10])
            legend('|\psi(R,t)|^2','V(R,t)')
            drawnow
        end
    end
end
%% ===================== 4 get all eigenstates ======================
% only needs to run once (and every time r changes) and it takes a bit time
% ================= evolution  ================
% all eigenstates
N_bs = 24; % bs = bound states
eigen = zeros(N_bs,Nr);
for N0 = 0:(N_bs -1)
    Ni = N0 +1;
    eigen_N0 = exp(-z/2).*z.^(lam -N0 -1/2) .*laguerreL(N0,2*lam-2*N0 -1,z);
    eigen(Ni,:) = eigen_N0/(sqrt(dr) *norm(eigen_N0));
end

%% ===================== 5 analysis =====================
% decompose into eigenstates
Nt_cycle_01 = size(psi,1);
dt_cycle_01 = dti_cycle_01/5 *dt;
t_cycle_01 = (1:Nt_cycle_01) *dt_cycle_01;

Ct = zeros(Nt_cycle_01,N_bs); % overlap coefficient at time t of psi(t) with the eigen states of V
Ct =  dr * (conj(eigen) * psi.').';
% Ct is a vectorized version of:
%for ti = 1:Nt_cycle_01
%    for Ni = 1:N_bs
%        Ct(ti,Ni) = sum(conj(eigen(Ni,:)).*psi(ti,:))*dr;
%    end
%end

% compute probabilities
Pt = zeros(Nt_cycle_01,N_bs +1);
Pt(:,1:(end-1)) = abs(Ct).^2;
Pt(:,end) = 1 -sum(Pt,2);
%% =============== plot E(t) and evolution of probabilities (Fig 1) ==========

if plot_flag
    figure('Position',[50 50 500 700]);
    tiledlayout(2,1,'Padding','tight')

    % probabilities in each eigenstate
    nexttile
    vec = [0,1,8,14,24]; % similar to fig 1a of the paper
    for vi = 1:length(vec)
        plot(t_cycle_01/cycle_01,Pt(:,vec(vi)+1),'LineWidth',2); hold on;
        leg_Str{vi}=['P_{' num2str(vec(vi)) '}'];
    end

    title('ES Population');
    xlabel('Time (cycles)');
    leg = [leg_Str(1:end-1), 'P_{diss}'];
    legend(leg,'NumColumns',round(vi/2));
    set(gca,'FontSize',16)

    nexttile
    t_c01 = t/cycle_01; % dimless
    plot(t_c01, U, t_c01, omega/omega_01,'LineWidth',2); hold on
    plot(t_c01, U.*cos(omega.*t))
    xlabel('Time (cycles)');
    ylim([-1.1,1.1]);
    title('Laser Pulse');
    legend('U(t)','\omega/\omega_{0,1}','E(t)','Location','southwest');
    set(gca,'FontSize',16)
end