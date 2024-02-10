% Reference:
% [1]  1990 Chelkowski, Efficient Molecular Dissociation by a Chirped 
%      Ultrashort Infrared Laser Pulse
%
%
%% ========================= 1 setup pulse ===============================
% user inputs 
clear
showplot = 1;

% ================= conversions and constants ====================

% conversions (atomic units to ...)
au2amu = 1/1822.888; % m_e to amu
au2ang = 0.52917721092; % a_0 to angstrom 
au2eV = 27.211385; % E_h to eV
au2as = 24.18884326505; % hbar/E_h to attosecond
au2Vpang = 51.4220652; % E_h/(e *a_0) to V/angstrom
au2D = 2.541746; % e *a_0 to Debye

% fine structure constant alpha = k_e *e^2/(hbar *c)
% Bohr_radius = hbar^2/(k_e *m_e *e^2) = hbar/(m_e*c*alpha)
% Hartree E_h = k_e^2 *m_e *e^4/hbar^2 = alpha^2 *m_e *c^2

% other conversion
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


% ================= Pulse ==========================
%I = 1e13; % W/cm^2
%E_M = sqrt(8*pi*k_e_SI/c_SI *I *1e4) *1e-10; % V/angstrom
%E_M = E_M/au2Vpang; % au = E_h/(e*a_0)


m_Ir = 192.217; % amu
m = m_Ir*m_Ir/(m_Ir +m_Ir); % amu
m = m/au2amu; % au, reduced mass

a = 1;%0.71741/a_0; % au = 1/a_0
D = 1;%6.125/au2eV; % au = E_h
%B = hbar *a/sqrt(2*m*D); % dimless
%omega_0 = 2*B*D/hbar; % au = E_h/hbar
%omega_01 = omega_0 *(1 -B); % au = E_h/hbar
%cycle_01 = 2*pi/omega_01; % au = hbar/E_h

N = 8; % dimless
alpha_0 = 2.5; % dimless
alpha_F = 6.25; % dimless
S = 1.5*pi; S_0 = 1.5*pi; % dimless
d_1 = 0.786/(a_0*au2D); % au = e
p_01 = 0.097/au2D; % au = e*a_0
 
% time steps
cycle_01=0.01;
dt = 0.001;% *cycle_01; % au = hbar/E_h
t_max = 10 *cycle_01; % au = hbar/E_h
t = 0:dt:t_max;

Nt = length(t);
t_max = max(t); % au = hbar/E_h
 
%% ======================= 2 initial state ============================
% ================= Morse Potential =================
r_0 = 2.9/au2ang;%1.7329*a_0; % au = a_0
Nr = 2^11;
r_max = 130; % au = a_0
r = linspace(0.1,r_max,Nr); % au = a_0
dr = r(2) -r(1);
x = r - r_0; % au = a_0
dx = dr;

% Morse potential at time t = 0
V = D*(1 - exp(-a*x)).^2; % E_h

 % ======================= ground state ===================
%lam = 1/B; % dimless
%z = 2*lam *exp(-a*x);

%n = 0;
psi_0 = exp(-(r-3/au2ang).^2/0.01); %z.^(lam -n -1/2) .*exp(-z/2) .*laguerreL(n,2*lam-2*n -1,z);
psi_0 = psi_0/(sqrt(dr)*norm(psi_0));
 
%% ================= 3 propagation (Split Step) ============
dti_cycle_01 = 0.01*round(cycle_01/(2*dt)); % index



dp = 2*pi*hbar/(Nr*dr); % au = hbar/a_0
p = ((0:Nr-1) -Nr/2) *dp; % au = hbar/a_0
T = p.^2/(2*m); % au = E_h
T = ifftshift(T); % au = E_h

dt_hbar = dt/hbar;
psi_1 = psi_0;

count = 0;
for ti = 1:Nt
    %V = D*(1 - exp(-a*x)).^2 - x *d_1*E_M*U(ti)*cos(omega(ti)*t(ti));
    
    psi_1 = exp(-1i *dt_hbar/2 *V) .*psi_1;
    psi_1 = ifft(exp(-1i *dt_hbar *T) .*fft(psi_1));
    psi_1 = exp(-1i *dt_hbar/2 *V) .*psi_1;
    
    psi_1 = psi_1/(sqrt(dr)*norm(psi_1)); % normalize
    
    if rem(ti,dti_cycle_01) == 0
        count = count +1;
        psi(count,:) = psi_1;
%         disp(ti);
        
        figure(89);
        plot(r, abs(psi_1).^2, r, V);
        title([num2str(ti*dt/cycle_01) ,' cycles']);
        ylim([0,10]); xlim([1/au2ang,7/au2ang])
        pause(0.02);
    end
end
 






%% ===================== 7 custom functions =====================
function addLine(pos, varargin)
% ADDLINE adds a line intersection given position.
%   pos: 1-by-2 row vector or #
%   
%
% Optional Inputs
%   XorY: stands for vertical ('X') or horizontal ('Y') line. By default,
%   == 'X'
%   Line Specs: see MATLAB plot for details
%
% Remarks
%   cannot change LineWidth for some reason
%
% Internal Func
%   addtick
%
%   last updated 02/05/2018
%

%% optional inputs
input = {'X', ''};
input_N = length(input);
for jj = 1:input_N
    if nargin >= 1+jj && ~isempty(varargin{jj})
        input{jj} = varargin{jj};
    end
end

[XorY, lineSpec] = cell2var(input);
%% main
Xm = get(gca, 'xlim'); % m stands for min and max
Ym = get(gca, 'ylim');


N = 50;

if length(pos) == 1
    if strcmp(XorY, 'Y')
        X = linspace(Xm(1), Xm(2), N);
        Y = pos*ones(1,N);
    elseif strcmp(XorY, 'X')
        X = pos*ones(1,N);
        Y = linspace(Ym(1), Ym(2), N);
    else
        error('invalid second input');
    end
    tick = pos;
elseif length(pos) == 2
    if strcmp(XorY, 'Y')
        X = linspace(Xm(1), pos(1), N);
        Y = pos(2)*ones(1,N);
        tick = pos(2);
    elseif strcmp(XorY, 'X')
        X = pos(1)*ones(1,N);
        Y = linspace(Ym(1), pos(2), N);
        tick = pos(1);
    else
        error('invalid second input');
    end
else
    error('invalid first input');
end
        
hold on;
plot(X, Y, lineSpec);
xlim(Xm); ylim(Ym); % in case the figure changes x,y limits
addTick(tick, XorY);

hold off;


end














% ================ custom functions ================
function addTick(tick1, varargin)
% ADDTICKS add ticks to a figure
%   tick: # or vector of tick values
%
% Optional Inputs
%   XorY: input 'X' or 'Y' to choose the axis on which the ticks is added. 
%   By default, == 'X'
%   ax: axes handle. By default, == gca (get current axes)
%   
% Custum Functions
%   cell2var
%
%   last updated 08/09/2017

%% optional inputs
input = {'X', gca};
input_N = length(input);
for jj = 1:input_N
    if nargin >= 1+jj && ~isempty(varargin{jj})
        input{jj} = varargin{jj};
    end
end

[XorY, ax] = cell2var(input);
%% main
XorYTick = [XorY, 'Tick']; 
tick0 = get(ax, XorYTick);
tick = unique([tick0, tick1]);
set(gca, XorYTick, tick);

end
function deleteTick(tick1, varargin)
%ADDTICKS add ticks to a figure
%   tick: # or vector of tick values
%
%Optional Inputs
%   XorY: input 'X' or 'Y' to choose the axis on which the ticks is deleted. 
%   By default, == 'X'
%   ax: axes handle. By default, == gca (get current axes)
%   
%Custum Functions
%   cell2var
%
%   last updated 08/09/2017

%% optional inputs
input = {'X', gca};
input_N = length(input);
for jj = 1:input_N
    if nargin >= 1+jj && ~isempty(varargin{jj})
        input{jj} = varargin{jj};
    end
end

[XorY, ax] = cell2var(input);
%% main
XorYTick = [XorY, 'Tick'];
N = length(tick1);
ticks = get(ax, XorYTick);
for kk = 1:N
    ticks(ticks == tick1(kk)) = [];
end
set(gca, XorYTick, ticks);

end












% ================ custom functions ================
function [varargout] = cell2var(cellA, varargin)
%CELL2VAR stands for cell to variable. It converts a cell of length N to N
%specified variables
%   cellA: 1D cell
%
%Optional Inputs
%   num: a # specifying which cell to convert to variable
%
%   last updated 06/05/2017
%   version 1.0
%
%   Author: Andrew Yuan
%   Copyright (c) 2017, All Rights Reserved
%

%% possible case
if nargin >= 2 && ~isempty(varargin{1})
    varargout{1} = cellA{varargin{1}};
    return;
end
%% main
if nargout <= length(cellA)
   for kk = 1:nargout
       varargout{kk} = cellA{kk};
   end
else
    error('number of outputs exceeds length of input cell');
end


end
