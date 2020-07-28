%% Problems 1 and 2
close all
%Initial variables and equations given in project document
x_00 = [0;00;100;100];
P_00 = [0 0 0 0; 0 0 0 0; 0 0 25 5; 0 0 5 25];
A = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1];
Bu_t = [0;0;0;-9.8];
Q =[20 0 5 0; 0 20 0 5; 5 0 10 0; 0 5 0 10];

%Initializes discrete time vector and sets up matrix to hold probability
%values for the cdf plot
figure
t = 1:1:30;
cdf_vect = zeros(1,30);
muy_vect = zeros(1,30);
mux_vect = zeros(1,30);
stdx_vect = zeros(1,30);
stdy_vect = zeros(1,30);
rho_vect = zeros(1,30);

%Same code as from Problem 1, but this iteration calculates the probability
%that the projectile has hit the ground by every time t <= 30 and puts that
%probability into a matrix
for i = 1:30
    x_00 = A*x_00 + Bu_t;
    P_00 = A*P_00*A.' + Q;
    prob_cdf = normcdf(0, x_00(2,1), sqrt(P_00(2,2)));
    cdf_vect(i) = prob_cdf;
    
    %keeps track of values needed for equations for 2d
    mux_vect(i) = x_00(1,1);
    muy_vect(i) = x_00(2,1);
    stdx_vect(i) = sqrt(P_00(2,2));
    stdy_vect(i) = sqrt(P_00(1,1));
    rho_vect(i) = (P_00(1,2)/(x_00(2,1)*x_00(1,1)));
    
    %Plots the vertical lines to data points
    plot([i i],[0 prob_cdf])
    hold on
end

%Plots the probabilities calculated above against their time stamps to form
%the discrete cdf of T_land
scatter(t,cdf_vect)
title('T_{land} CDF')
xlabel('Time')
ylabel('CDF')
figure

%This block of code finds the height differences between the cdf points,
%which are the probabilities of the points (pmf) since this is a discrete system.
pmf_vect = [cdf_vect(1)];
for j = 2:30
    prob_pmf = cdf_vect(j) - cdf_vect(j-1);
    pmf_vect = [pmf_vect; prob_pmf];
    plot([j j],[0 prob_pmf])
    hold on
end

%Plots the pmf versus time
scatter(t, pmf_vect)
title('T_{land} PMF')
xlabel('Time')
ylabel('Probability')
pmf_vect_t = [pmf_vect.';t];

%Initializes vectors for containing the summation values
x = linspace(1000,3500);
values = 0;
final = 0;
new_mu = zeros(1,30);
new_std = zeros(1,30);
figure
sum_mean = 0;
sum_Ex2 = 0;

%Looking at all discrete t from 1 to 30, every 1 second
for a = 1:30
    %Calculates the various Marginal pdfs, evaluates them over the specified x
    %interval and scales them with their respective probabilities, plots the 
    %individuals, then sums their values into the final variable
    new_mu(a) = mux_vect(a) + (stdx_vect(a)/stdy_vect(a)*rho_vect(a)*(-muy_vect(a)));
    new_std(a) = sqrt((1-rho_vect(a)^2)*stdx_vect(a)^2);
    values = (pmf_vect(a)*normpdf(x,new_mu(a),new_std(a)));
    final = final+values;
    title('Various Marginal PDFs')
    xlabel('x^1 Distance')
    ylabel('PDF')
    plot(x, values)
    hold on
    
    %Finds the expected values of the marginal pdfs separately, then adds
    %them together to get the expected value of the mixture
    fun_Ex = @(y) (y.*pmf_vect(a).*normpdf(y,new_mu(a),new_std(a)));
    individual_mean = integral(fun_Ex,1000,3500);
    sum_mean = sum_mean + individual_mean;
    
    %Finds E[x^2] for each marginal pdf, then adds them together to get the
    %mixture value
    fun_Ex2 = @(y) (y.^2.*pmf_vect(a).*normpdf(y,new_mu(a),new_std(a)));
    individual_Ex2 = integral(fun_Ex2,1000,3500);
    sum_Ex2 = sum_Ex2 + individual_Ex2;
end

%Prints the mixture mean, calculates its variance and then prints it as
%well
sum_mean
sum_variance = sum_Ex2 - (sum_mean)^2

%Plots the summation of the marginals scaled by probability
figure
plot(x, final)
title('f_X_{land}(x)')
xlabel('x^1 Distance')
ylabel('PDF')

%% Problem 3 base vs radar
clear all
close all

%Initial variables and equations given in project document
x_pre = [0;00;100;100];
P_pre = [0 0 0 0; 0 0 0 0; 0 0 25 5; 0 0 5 25];
A = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1];
Bu_t = [0;0;0;-9.8];
Q =[20 0 5 0; 0 20 0 5; 5 0 10 0; 0 5 0 10];
w_t_mu=[0,0,0,0];
R = [9 3 ; 3 9];
r_t_mu = [0 0];
C = [1 0 0 0;0 1 0 0];

%Initializes discrete time vector and sets up vectors to hold position
%values
t = 1:1:30;
x_position_vect = [];
x2_position_vect = [];
x_radar_vect = [];
y_radar_vect = [];

%Generates an initial condition from the bivariate Gaussian, transposes to
%have correct format
x_t = mvnrnd(x_pre,P_pre,1);
x_t = x_t.';
    
%For 30 timesteps, generates a random w_t value and uses that, with the
%given equations, to predict the position at the next timestep.  Loop is
%broken if projectile hits the ground before 30
for i = 1:30
    %This section calculates the radar measurements y 
    r_t = mvnrnd(r_t_mu,R,1).';
    y_t = C*x_t + r_t;
    
    w_t = mvnrnd(w_t_mu,Q,1);
    x_t = A*x_t + Bu_t + w_t;
    
    if x_t(2) <= 0
        break
    end
 
    x_position_vect = [x_position_vect x_t(1)];
    x2_position_vect = [x2_position_vect x_t(2)];
    
    x_radar_vect(i) = y_t(1);
    y_radar_vect(i) = y_t(2);
         
end
   
%Plots the radar observations versus the ground truth measurements  
figure()
scatter(x_position_vect,x2_position_vect)
hold on
xlabel('x^1 Distance')
ylabel('x^2 Distance')
title('One Realization of Trajectory with Radar')
scatter(x_radar_vect,y_radar_vect)
legend('Ground Truth','Radar Observations')
hold off

%% Problem 3 Kalman filter and graph comparison

close all
%Initial variables and equations given in project document
x_pre = [0;00;100;100];
P_pre = [0 0 0 0; 0 0 0 0; 0 0 25 5; 0 0 5 25];
A = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1];
Bu_t = [0;0;0;-9.8];
Q =[20 0 5 0; 0 20 0 5; 5 0 10 0; 0 5 0 10];
w_t_mu=[0,0,0,0];
R = [9 3 ; 3 9];
r_t_mu = [0 0];
C = [1 0 0 0;0 1 0 0];

%Generates an initial condition from the bivariate Gaussian, transposes to
%have correct format
x_t = mvnrnd(x_pre,P_pre,1);
x_t = x_t.';

%Running kalman filter for 10 seconds
for i = 1:10

    %This section calculates the radar measurements y 
    r_t = mvnrnd(r_t_mu,R,1).';
    y_t = C*x_t + r_t;

%Updates the state x
    w_t = mvnrnd(w_t_mu,Q,1).';
    x_t = A*x_t + Bu_t + w_t;

%Kalman filtering subscript pre for prediction, up for update 
    x_pre = A*x_pre + Bu_t;
    P_pre = A*P_pre*A.' + Q;
    y_up = y_t - C*x_pre;
    S_up = C*P_pre*C.' + R;
    K_up = P_pre*C.'*S_up^-1;
    x_up = x_pre + K_up*y_up;
    P_up = (eye(4)-K_up*C)*P_pre;
end

x_t_kalman = x_up
P_t_kalman = P_up

%Initial matrices
x_up = [0;00;100;100];
P_up = [0 0 0 0; 0 0 0 0; 0 0 25 5; 0 0 5 25];
A = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1];
Bu_t = [0;0;0;-9.8];
Q =[20 0 5 0; 0 20 0 5; 5 0 10 0; 0 5 0 10];

%Updates value for five time steps of one second
for i = 1:10
        x_up = A*x_up + Bu_t;
        P_up = A*P_up*A.' + Q;
end

%Outputs final results and corrects variable name
x_100 = x_up
P_100 = P_up

figure
hold on
confidence_ellipse(x_100,P_100, .9, 'r')
confidence_ellipse(x_t_kalman,P_t_kalman, .9, 'b')
title('90% CI at t = 10 With and Without Filtering')
xlabel('x^1 Distance')
ylabel('x^2 Distance')
legend('No Filter','Kalman Filter')
hold off
%% Problem 3 kalman up to t = 10 then prediction without kalman
close all
%Initial variables and equations given in project document
x_pre = [0;00;100;100];
P_pre = [0 0 0 0; 0 0 0 0; 0 0 25 5; 0 0 5 25];
A = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1];
Bu_t = [0;0;0;-9.8];
Q =[20 0 5 0; 0 20 0 5; 5 0 10 0; 0 5 0 10];
w_t_mu=[0,0,0,0];
R = [9 3 ; 3 9];
r_t_mu = [0 0];
C = [1 0 0 0;0 1 0 0];

%Generates an initial condition from the bivariate Gaussian, transposes to
%have correct format
x_t = mvnrnd(x_pre,P_pre,1);
x_t = x_t.';

%Sets up time vector and a vector to store points for the cdf
figure
t = 1:1:30;
cdf_vect = zeros(1,30);
cdf_vect = zeros(1,30);
muy_vect = zeros(1,30);
mux_vect = zeros(1,30);
stdx_vect = zeros(1,30);
stdy_vect = zeros(1,30);
rho_vect = zeros(1,30);
x_position_vect = [];
x2_position_vect = [];
x_radar_vect = [];
y_radar_vect = [];

%Running kalman filter for 10 seconds
for i = 1:10

    %This section calculates the radar measurements y 
    r_t = mvnrnd(r_t_mu,R,1).';
    y_t = C*x_t + r_t;

%Updates the state x
    w_t = mvnrnd(w_t_mu,Q,1).';
    x_t = A*x_t + Bu_t + w_t;

%Kalman filtering subscript pre for prediction, up for update 
    x_pre = A*x_pre + Bu_t;
    P_pre = A*P_pre*A.' + Q;
    y_up = y_t - C*x_pre;
    S_up = C*P_pre*C.' + R;
    K_up = P_pre*C.'*S_up^-1;
    x_up = x_pre + K_up*y_up;
    P_up = (eye(4)-K_up*C)*P_pre;
    
    %Calculates the probability that the projectile has hit the ground by
    %the given time and stores it in the vector
    prob_cdf = normcdf(0, x_up(2,1), sqrt(P_up(2,2)));
    cdf_vect(i) = prob_cdf;
    
    %keeps track of values needed for equations for position calculation
    mux_vect(i) = x_up(1,1);
    muy_vect(i) = x_up(2,1);
    stdx_vect(i) = sqrt(P_up(2,2));
    stdy_vect(i) = sqrt(P_up(1,1));
    rho_vect(i) = (P_up(1,2)/(x_up(2,1)*x_up(1,1)));
    
    %Plots the vertical lines to data points
    plot([i i],[0 prob_cdf])
    hold on
end

kalman_mean = x_up
kalman_covariance = P_up

%Same code as from problem 2, but this only starts after the radar has
%stopped working and there are no measurements
for i = 11:30
    x_up = A*x_up + Bu_t;
    P_up = A*P_up*A.' + Q;
    prob_cdf = normcdf(0, x_up(2,1), sqrt(P_up(2,2)));
    cdf_vect(i) = prob_cdf;
    
    %keeps track of values needed for equations for position calculation
    mux_vect(i) = x_up(1,1);
    muy_vect(i) = x_up(2,1);
    stdx_vect(i) = sqrt(P_up(2,2));
    stdy_vect(i) = sqrt(P_up(1,1));
    rho_vect(i) = (P_up(1,2)/(x_up(2,1)*x_up(1,1)));
    
    
    %Plots the vertical lines to data points
    plot([i i],[0 prob_cdf])
    hold on
end

%Plots the probabilities calculated above against their time stamps to form
%the discrete cdf of T_land
scatter(t,cdf_vect)
title('T_{land} CDF')
xlabel('Time')
ylabel('CDF')
figure


%This block of code finds the height differences between the cdf points,
%which are the probabilities of the points (pmf) since this is a discrete system.
pmf_vect = [cdf_vect(1)];
for j = 2:30
    prob_pmf = cdf_vect(j) - cdf_vect(j-1);
    pmf_vect = [pmf_vect; prob_pmf];
    plot([j j],[0 prob_pmf])
    hold on
end

%Plots the pmf versus time
scatter(t, pmf_vect)
title('T_{land} PMF')
xlabel('Time')
ylabel('Probability')
pmf_vect_t = [pmf_vect.';t];
ylim([0,0.5]);

%Initializes vectors for containing the summation values
x = linspace(1000,3500);
values = 0;
final = 0;
new_mu = zeros(1,30);
new_std = zeros(1,30);
figure
sum_mean = 0;
sum_Ex2 = 0;

%Looking at all discrete t from 1 to 30, every 1 second
for a = 1:30
    %Calculates the various Marginal pdfs, evaluates them over the specified x
    %interval and scales them with their respective probabilities, plots the 
    %individuals, then sums their values into the final variable
    new_mu(a) = mux_vect(a) + (stdx_vect(a)/stdy_vect(a)*rho_vect(a)*(-muy_vect(a)));
    new_std(a) = sqrt((1-rho_vect(a)^2)*stdx_vect(a)^2);
    values = (pmf_vect(a)*normpdf(x,new_mu(a),new_std(a)));
    final = final+values;
    title('Various Marginal PDFs')
    xlabel('x^1 Distance')
    ylabel('PDF')
    plot(x, values)
    hold on
    
    %Finds the expected values of the marginal pdfs separately, then adds
    %them together to get the expected value of the mixture
    fun_Ex = @(y) (y.*pmf_vect(a).*normpdf(y,new_mu(a),new_std(a)));
    individual_mean = integral(fun_Ex,1000,3500);
    sum_mean = sum_mean + individual_mean;
    
    %Finds E[x^2] for each marginal pdf, then adds them together to get the
    %mixture value
    fun_Ex2 = @(y) (y.^2.*pmf_vect(a).*normpdf(y,new_mu(a),new_std(a)));
    individual_Ex2 = integral(fun_Ex2,1000,3500);
    sum_Ex2 = sum_Ex2 + individual_Ex2;
end

%Prints the mixture mean, calculates its variance and then prints it as
%well
sum_mean
sum_variance = sum_Ex2 - (sum_mean)^2

%Plots the summation of the marginals scaled by probability
figure
plot(x, final)
title('f_X_{land}(x)')
xlabel('x^1 Distance')
ylabel('PDF')
%% Plotting possible realizations after radar stops
close all
%Initial variables and equations given in project document, with x10/10 and
%p10/10 calculated from Kalman filtering in part c
x_10 = [909.7296,582.0647,88.6520,5.0765];
P_10 = [8.9857 2.9922 1.1185 0.3557; 2.9922 8.9857 0.3557 1.1185; 1.1185 0.3557 31.3630 0.3569; 0.3557 1.1185 0.3569 31.3630];
x_00 = [0;00;100;100];
P_00 = [0 0 0 0; 0 0 0 0; 0 0 25 5; 0 0 5 25];
A = [1 0 1 0; 0 1 0 1; 0 0 1 0; 0 0 0 1];
Bu_t = [0;0;0;-9.8];
Q =[20 0 5 0; 0 20 0 5; 5 0 10 0; 0 5 0 10];
w_t_mu=[0,0,0,0];

%Initializes discrete time vector and sets up matrix to hold probability
%values for the cdf plot
t = 1:1:30;
x_position_vect = [];
y_position_vect = [];
figure
P_t=P_10;

%Cell arrays to hold data for average position and time of landing calculation
x_pos={};
ts={};

%Plots 1000 realizations of possible trajectories after the radar has
%stopped working at t = 10
for n = 1:3000
    %Generates an initial condition from the bivariate Gaussian, transposes to
    %have correct format
    x_t = mvnrnd(x_10,P_10,1);
    x_t = x_t.';
    
    %For 15 timesteps past the radar turn off, generates a random w_t value
    %and uses that, with the given equations, to predict the position at 
    %the next timestep. Loop is broken if projectile hits the ground before 30
    for i = 1:15
        w_t = mvnrnd(w_t_mu,Q,1).';
        x_t = A*x_t + Bu_t + w_t;
        P_t = A*P_t*A.' + Q;
    
        %When projectile hits the ground, position and time are recorded
        if x_t(2,1) <= 0
            x_pos=[x_pos x_t(1)];
            ts=[ts i];
            break
        end
 
        x_position_vect = [x_position_vect x_t(1)];
        y_position_vect = [y_position_vect x_t(2)];
         
    end
   
    
    %Plots the specific realization and clears vector to hold data for the
    %next
    scatter(x_position_vect,y_position_vect)
    hold on
    x_position_vect = [];
    y_position_vect = [];

    
end

%Plots the various trajectories and calculates the average landing distance
%and time based on the data from the trials
xlabel('x^1 Distance')
ylabel('x^2 Distance')
title('Final Trajectory')
x_land=sum([x_pos{:}])/length(x_pos)
t_land=sum([ts{:}])/length(ts) + 10

