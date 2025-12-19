function GKM = gkm_numeric(param,timeInfo)

%--------------------------------------------------------------------------
% Input parameters:
%   alpha           % inversion efficacy
%   T1              % tissue T1
%   Delta_t         % transit delay
%   tau             % perfusion duration
%   m               % magnetization
%   f               % CBF ml of perfusive substance per ml tissue per sec
%   lambda          % partition coefficient water/tissue
%--------------------------------------------------------------------------
% vmalis@ucsd.edu
% Vadim Malis

T1=param(1);
Delta_t=param(2);
tau=param(3);
f=param(4);

lambda=1;
alpha=1;

ns_time=5;
tag=20;
dt=0.1;
t_total=5000;

% Bloch
[mz,~]=blochMz(T1,1,dt,t_total,[ns_time,tag]);
control=ones(size(mz,1),1);
    
% GKM
m = control - mz(:,1,end-1);
%m = m(ns_time/dt:(t_total-max(tag)+ns_time)/dt-1);

%t=0:dt:(t_total/dt-tag/dt-dt)*dt;
t=0:dt:t_total-dt;

% Define the piecewise function c(t) for the pulsed state
c = zeros(size(t));
c(t >= Delta_t & t < tau + Delta_t) = alpha * exp(-t(t >= Delta_t & t < tau + Delta_t)/T1);

% Define r(t) and m(t)
r  = exp(-f*t/lambda);

% r(t) and m(t)
rm = r.*m'; % Multiplying

% Convolve c(t) with rm_conv
Delta_M_conv = f * conv(c, rm, 'full')*dt; % Multiply by dt to approximate the integral

% Truncate the convolution result to match the time vector length
GKM = Delta_M_conv(1:length(t));
GKM(1:ceil(Delta_t/dt)) = 0;
GKM = interp1(t, GKM, timeInfo, 'linear');

end