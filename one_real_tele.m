%Helene Levy
%MAE 271B Project
%Telegraphing
clc; clear; close all;

%% Given Parameters and Statistics
% time parameters
dt = 0.001;
tf = 10; % sec
tau = 2; % sec

t = 0:dt:tf;

% lateral position
m_y = 0;
var_y = 0;

% lateral velocity 
m_v = 0;
var_v = 200^2; % (ft/sec)^2
corr_yv = 0;

% target acceleration
m_at = 0;
var_at = 100^2; % (ft/sec^2)^2

% target acceleration
m_aTbar = 0;
aT = 100; %ft/sec^2
RaT = @(t,s) aT^2*exp(-2*lambda*abs(t-s));

%relative velocity
R1 = 15*10^(-6); %rad^2/sec
R2 = 1.67*10^(-3); %rad^2/sec^3
Vc = 300; %ft/sec

% process noise spectral density
W = [ 0 0 0; 0 0 0; 0 0 var_at];

%fading and scintillation noise
m_n = 0;
V = @(t) R1 + R2/(tf-t)^2 ; 

% initial covariance
P0 = [var_y 0 0; 0 var_v 0; 0 0 var_at];

%% Telegraph Process Generation
%poisson variable param
lambda = 0.25; %/sec

a =[];

%initial aT 
A = rand; 
if A <= 0.5
    a(1) = aT;
else
    a(1) = -aT;
end

%generating random switching times for telegraph signal 
time = [];
time(1) = 0; i = 1;
t_np1 = 0;
%only want times in our simulation <10
while t_np1 < 10 
    %uniform random variable
    U = rand;
    %next random switching time
    t_np1 = time(i) - 1/lambda*log(U);
    if t_np1 < 10
        time(i+1) = t_np1;
        a(i+1) = -a(i);
    end
    i = i+1;
end

%rounding the random time to accuracy of dt
time = round(time,-log10(dt));
t = round(t,-log10(dt));
[~,loc] = ismember(time,t);

%converting a to timescale 
at = zeros(1,length(t));
for j = 2:length(loc)
    at(loc(j-1):loc(j)) = a(j-1);
end
at(loc(end):end) = a(end);

%visualizing telegraph signal
stairs(t,at)
ylim([-110,110])

%% Kalman Filter
% State space matrices
%two state telegraph
Ft = [0 1; 0 0];
Gt = [0; -1];

%original three state
F = [0 1 0; 0 0 -1; 0 0 -1/tau]; 
B = [0; 1; 0];
G = [0; 0; 1];
H = @(t) [1/(Vc*(tf-t)) 0 0];
Hbar = [1 0 0];

%initial states
y_0 = m_y + sqrt(var_y)*randn(1,1); %(0)
v_0 = m_v + sqrt(var_v)*randn(1,1);

% noise generation
n = @(t) m_n + sqrt(V(t)/dt)*randn(1);

X0 = [y_0 v_0 at(1)]';

%preallocation matrices
Xhat = zeros(3,length(t));
X = zeros(3,length(t));
P = zeros(3,3,length(t));
z = zeros(1,length(t)-1);
r = zeros(1,length(t)-1);

K1 = zeros(1,length(t));
K2 = zeros(1,length(t));
K3 = zeros(1,length(t));

X(:,1) = X0;
P(:,:,1) = P0;
K1(1) = P(1,1,1)/(Vc*R1*(tf-t(1)) + Vc*R2/(tf-t(1)));
K2(1) = P(1,2,1)/(Vc*R1*(tf-t(1)) + Vc*R2/(tf-t(1)));
K3(1) = P(1,3,1)/(Vc*R1*(tf-t(1)) + Vc*R2/(tf-t(1)));
a_p = 0;

for i = 1: length(t)-1
    %true states
    dX = (Ft*X(1:2,i)+Gt*at(i))*dt;
    X(1:2,i+1) = X(1:2,i) + dX;
    X(3,i+1) = at(i+1);
    
    %measurement
    z(i) = H(t(i))*X(:,i)+n(t(i));
    
    %covariance matrix
    P_i = P(:,:,i);
    Pdot = F*P_i + P_i*transpose(F) - 1/(Vc^2*R1*(tf-t(i))^2+Vc^2*R2)*P_i...
     *Hbar'*Hbar*P_i+W;
    P(:,:,i+1) = Pdot*dt + P(:,:,i);
 
    %calculating Kalman gains
    K1(i+1) = P(1,1,i+1)/(Vc*R1*(tf-t(i+1)) + Vc*R2/(tf-t(i+1)));
    K2(i+1) = P(1,2,i+1)/(Vc*R1*(tf-t(i+1)) + Vc*R2/(tf-t(i+1)));
    K3(i+1) = P(1,3,i+1)/(Vc*R1*(tf-t(i+1)) + Vc*R2/(tf-t(i+1)));
    
    %estimates
    K = [K1(i) K2(i) K3(i)]';
    dXhat = F*Xhat(:,i)*dt + K*(z(i)-H(t(i))*Xhat(:,i))*dt;
    Xhat(:,i+1) = Xhat(:,i) + dXhat;
end
Ehat = X - Xhat;
sig= sqrt([squeeze(P(1,1,:)), squeeze(P(2,2,:)), squeeze(P(3,3,:))]);

%% Plotting Results
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

