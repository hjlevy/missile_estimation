%Helene Levy
%MAE 271B Project
%Monte Carlo Simulation
clc; clear; close all;

%time step choices
dt = 0.001;
tf = 10;
t = 0:dt:tf;

%preallocations
e_sum = zeros(3,length(t));
e_ave = zeros(3,length(t));
r_sum = zeros(1,1);

P_ave = zeros(3,3,length(t));
P_sum = zeros(3,3,length(t));

%iterating through N times
N = 10000;
for j = 1:N
    [X_hat,e_hat,P,r] = tele_kalman_filt(dt);
    
    %ensemble average of error 
    e_sum = e_sum + e_hat;
    e_ave = 1/j*e_sum;
    
    for k = 1:length(e_hat)
        %actual error variance calculation
        P_sum(:,:,k) = P_sum(:,:,k)+(e_hat(:,k)-e_ave(:,k))...
                       *(e_hat(:,k)-e_ave(:,k))';
        P_ave(:,:,k) = 1/(j-1)*P_sum(:,:,k);
        
    end
    
    %ensemble average for correlation of residuals at end time
    r_sum = r_sum + r(end/2)*(r(end/2-1))';
    r_ave = 1/j*r_sum;
end

%showing residuals uncorrelated
disp('ensemble average for correlation of the residuals:');
disp(r_ave);

%apriori standard deviation and actual standard deviation
std_true = sqrt([squeeze(P(1,1,:)), squeeze(P(2,2,:)), squeeze(P(3,3,:))]);
std_ave = sqrt([squeeze(P_ave(1,1,:)), squeeze(P_ave(2,2,:)),...
                squeeze(P_ave(3,3,:))]);

%plotting std comparison
figure; 
plot(tf-t,std_true(:,1)); hold on;
plot(tf-t,std_ave(:,1));
set ( gca, 'xdir', 'reverse' );
xlabel('time to go (sec)');
ylabel('Standard deviation of the position error');
title('Position: Evolution of Actual Error RMS vs. Apriori')
legend('a priori','telegraph actual')

figure;
plot(tf-t,std_true(:,2)); hold on;
plot(tf-t,std_ave(:,2));
set ( gca, 'xdir', 'reverse' );
xlabel('time to go (sec)');
ylabel('Standard deviation of the velocity error');
title('Velocity: Evolution of Actual Error RMS vs. Apriori')
legend('a priori','telegraph actual')

figure;
plot(tf-t,std_true(:,3)); hold on;
plot(tf-t,std_ave(:,3));
set ( gca, 'xdir', 'reverse' );
xlabel('time to go (sec)');
ylabel('Standard deviation of the acceleration error');
title('Acceleration: Evolution of Actual Error RMS vs. Apriori')
legend('a priori','telegraph actual')

