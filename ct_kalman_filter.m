function [Xhat,Ehat,P,r] = ct_kalman_filter(dt)
% function for running continuous time kalman filter

% Given parameters and statistics

% time parameters
tf = 10; % sec
tau = 2; % sec

R1 = 15*10^(-6); %rad^2/sec
R2 = 1.67*10^(-3); %rad^2/sec^3
Vc = 300; %ft/sec

% target acceleration
m_at = 0;
var_at = 100^2; % (ft/sec^2)^2

% lateral position
m_y = 0;
var_y = 0;

% lateral velocity 
m_v = 0;
var_v = 200^2; % (ft/sec)^2

%fading and scintillation noise
m_n = 0;
V = @(t) R1 + R2/(tf-t)^2 ; 

% process noise spectral density
W = [ 0 0 0; 0 0 0; 0 0 var_at];

% initial covariance
P0 = [var_y 0 0; 0 var_v 0; 0 0 var_at];

% State Space Matrices
% xdot = Fc + Bap + Gw_at
% x = [y v at]'
F = [0 1 0; 0 0 -1; 0 0 -1/tau]; 
B = [0; 1; 0];
G = [0; 0; 1];

H = @(t) [1/(Vc*(tf - t)) 0 0];
Hbar = [1 0 0];

%K and P calculations using ode45
tspan = 0:dt:tf;
[t,P] = ode45(@Pdot, tspan , P0(:));
P = reshape(P.',3,3,[]);

%calculating Kalman gains
K1 = squeeze(P(1,1,:))./(Vc*R1.*(tf-t) + Vc*R2./(tf-t));
K2 = squeeze(P(1,2,:))./(Vc*R1.*(tf-t) + Vc*R2./(tf-t));
K3 = squeeze(P(1,3,:))./(Vc*R1.*(tf-t) + Vc*R2./(tf-t));

%initial states
y_0 = m_y + sqrt(var_y)*randn(1,1); %(0)
v_0 = m_v + sqrt(var_v)*randn(1,1);
at_0 = m_at + sqrt(var_at)*randn(1,1);

% noise generation
w_at = m_at + sqrt(var_at/dt)*randn(1,length(t)-1);
n = @(t) m_n + sqrt(V(t)/dt)*randn(1);

%preallocating matrices
Xhat = zeros(3,length(t));
X = zeros(3,length(t));
z = zeros(1,length(t)-1);
r = zeros(1,length(t)-1);

X0 = [y_0 v_0 at_0]';

%initial conditions of state
X(:,1) = X0;
a_p = 0;
for i = 1: length(t)-1
    %true states
    dX = (F*X(:,i)+B*a_p+G*w_at(i))*dt;
    X(:,i+1) = X(:,i) + dX;
    
    %measurement
    z(i) = H(t(i))*X(:,i)+n(t(i));
    
    %residuals
    r(i) = z(i)-H(t(i))*Xhat(:,i);
    
    %estimates
    K = [K1(i) K2(i) K3(i)]';
    dXhat = F*Xhat(:,i)*dt + K*(z(i)-H(t(i))*Xhat(:,i))*dt;
    Xhat(:,i+1) = Xhat(:,i) + dXhat;
    
end
Ehat = X - Xhat;
end