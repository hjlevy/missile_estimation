% MAE 271B Project
% Helene Levy
clc; clear; close all;

%% Given parameters and statistics

% time parameters
tf = 10; % sec
tau = 2; % sec

R1 = 15*10^(-6); %rad^2/sec
R2 = 1.67*10^(-3); %rad^2/sec^3
Vc = 300; %ft/sec

% target acceleration
m_at = 0;
var_at = 100^2; % (ft/sec^2)^2
corr_at_as = @(t,s) var_at*exp(-(t-s)/tau);

% lateral position
m_y = 0;
var_y = 0;

% lateral velocity 
m_v = 0;
var_v = 200^2; % (ft/sec)^2
corr_yv = 0;

%fading and scintillation noise
m_n = 0;
V = @(t) R1 + R2/(tf-t)^2 ; 
corr_nt_ntau = @(t) (R1 + R2/(tf-t)^2)*dirac(t-tau);

% process noise spectral density
W = [ 0 0 0; 0 0 0; 0 0 var_at];

% initial covariance
P0 = [var_y 0 0; 0 var_v 0; 0 0 var_at];


%% State Space Matrices
% xdot = Fc + Bap + Gw_at
% x = [y v at]'
F = [0 1 0; 0 0 -1; 0 0 -1/tau]; 
B = [0; 1; 0];
G = [0; 0; 1];

H = @(t) [1/(Vc*(tf - t)) 0 0];
Hbar = [1 0 0];

%% P, K calculation and Plotting Figures 9.6 and 9.7
dt = 0.0001;
tspan = 0:dt:tf;
[t,P] = ode45(@Pdot, tspan , P0(:));
P = reshape(P.',3,3,[]);

%calculating Kalman gains
K1 = squeeze(P(1,1,:))./(Vc*R1.*(tf-t) + Vc*R2./(tf-t));
K2 = squeeze(P(1,2,:))./(Vc*R1.*(tf-t) + Vc*R2./(tf-t));
K3 = squeeze(P(1,3,:))./(Vc*R1.*(tf-t) + Vc*R2./(tf-t));

%plotting Kalman gains vs. time
figure;
plot(t,K1,'b-'); hold on;
plot(t,K2,'r--'); hold on;
plot(t,K3,'m-.'); 

legend('K1','K2','K3');
xlabel('time (sec)');
ylabel('Kalman Filter Gain');
title('Filter Gain History');

figure;
rms_y = sqrt(squeeze(P(1,1,:)));
rms_v = sqrt(squeeze(P(2,2,:)));
rms_at = sqrt(squeeze(P(3,3,:)));
plot(t,rms_y,'b-'); hold on;
plot(t,rms_v,'r--'); hold on;
plot(t,rms_at,'m-.');

legend({'position (ft)','velocity (ft/sec)','acceleration (ft/sec^2)'});
xlabel('time (sec)');
ylabel('Standard deviation of the state error');
title('Evolution of Estimation Error RMS');

%% Kalman Filter
%initial states
y_0 = m_y + sqrt(var_y)*randn(1,1); %(0)
v_0 = m_v + sqrt(var_v)*randn(1,1);
at_0 = m_at + sqrt(var_at)*randn(1,1);

% noise generation
w_at = m_at + sqrt(var_at/dt)*randn(1,length(t)-1);
n = @(t) m_n + sqrt(V(t)/dt)*randn(1);

X0 = [y_0 v_0 at_0]';

Xhat = zeros(3,length(t));

X = zeros(3,length(t));
X(:,1) = X0;
a_p = 0;

for i = 1: length(t)-1
    %true states
    dX = (F*X(:,i)+B*a_p+G*w_at(i))*dt;
    X(:,i+1) = X(:,i) + dX;
    
    %measurement
    z = H(t(i))*X(:,i)+n(t(i));
   
    %estimates
    K = [K1(i) K2(i) K3(i)]';
    dXhat = F*Xhat(:,i)*dt + K*(z-H(t(i))*Xhat(:,i))*dt;
    Xhat(:,i+1) = Xhat(:,i) + dXhat;
end
Ehat = X - Xhat;
sig= sqrt([squeeze(P(1,1,:)), squeeze(P(2,2,:)), squeeze(P(3,3,:))]);

%plotting states
figure;
subplot(311)
plot(t,X(1,:)); hold on;
plot(t,Xhat(1,:))
legend('true', 'estimate');
title('Position');
xlabel('time (s)');

subplot(312)
plot(t,X(2,:)); hold on;
plot(t,Xhat(2,:))
legend('true', 'estimate');
title('Velocity');
xlabel('time (s)');

subplot(313)
plot(t,X(3,:)); hold on;
plot(t,Xhat(3,:))
legend('true', 'estimate');
title('Target Acceleration');
xlabel('time (s)');

%plotting state errors

figure;
subplot(311)
stairs(t,Ehat(1,:));hold on;
plot(t,sig(:,1),'r--',t,-sig(:,1),'r--'); 
title('Position Error');
xlabel('time (s)');
ylabel('position error');

subplot(312)
stairs(t,Ehat(2,:)); hold on;
plot(t,sig(:,2),'r--',t,-sig(:,2),'r--'); 
title('Velocity Error');
xlabel('time (s)');
ylabel('velocity error');

subplot(313)
stairs(t,Ehat(3,:)); hold on;
plot(t,sig(:,3),'r--',t,-sig(:,3),'r--'); 
title('Target Acceleration Error');
xlabel('time (s)');
ylabel('acceleration error');


