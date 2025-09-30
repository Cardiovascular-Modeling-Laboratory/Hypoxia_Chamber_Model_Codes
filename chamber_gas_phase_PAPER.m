%function [TR_10, TR_60, TR_90, BL_30, BL_60, BL_90] = chamber_gas_phase_PAPER(T1,td,Inlet_l2,Inlet_l3,cp,tf,O2_f)

%% Temperature data
% Find time on temperature curve using Well 1
if T1>=20.3 && T1<27.2
    t_m1 = (1800/6.9)*T1 - 36540/6.9;
elseif T1>=27.2 && T1<28.8
    t_m1 = (2700/1.6)*T1 - 70560/1.6;
elseif T1>=28.8 && T1<29.3
    t_m1 = (2700/0.5)*T1 - 75510/0.5;
else
    t_m1 = 7200;
end

% Determine gas temperatures above remaining wells and for Layers 1 and 2
% Well 2
if t_m1>=0 && t_m1<1800
    T2 = (7.6/1800)*t_m1 + 36720/1800;
elseif t_m1>=1800 && t_m1<4500
    T2 = (1.3/2700)*t_m1 + 73260/2700;
elseif t_m1>=4500 && t_m1<7200
    T2 = (0.1/2700)*t_m1 + 78660/2700;
else
    T2 = 29.4;
end

% Well 3
if t_m1>=0 && t_m1<1800
    T3 = (8.6/1800)*t_m1 + 36900/1800;
elseif t_m1>=1800 && t_m1<4500
    T3 = (1.6/2700)*t_m1 + 75690/2700;
elseif t_m1>=4500 && t_m1<7200
    T3 = -(0.1/2700)*t_m1 + 83340/2700; 
else
    T3 = 30.6;
end

% Well 4
if t_m1>=0 && t_m1<1800
    T4 = (9.3/1800)*t_m1 + 36720/1800;
elseif t_m1>=1800 && t_m1<4500
    T4 = (1.1/2700)*t_m1 + 78210/2700;
elseif t_m1>=4500 && t_m1<7200
    T4 = (0.4/2700)*t_m1 + 81360/2700;
else
    T4 = 31.2;
end

% Well 5
if t_m1>=0 && t_m1<1800
    T5 = (8.5/1800)*t_m1 + 20.4;
elseif t_m1>=1800 && t_m1<4500
    T5 = (1.1/2700)*t_m1 + 76050/2700;
elseif t_m1>=4500 && t_m1<7200
    T5 = (0.5/2700)*t_m1 + 78750/2700;
else
    T5 = 30.5;
end

% Well 6
if t_m1>=0 && t_m1<1800
    T6 = (8.6/1800)*t_m1 + 20.5;
elseif t_m1>=1800 && t_m1<4500
    T6 = (1.1/2700)*t_m1 + 76590/2700;
elseif t_m1>=4500 && t_m1<7200
    T6 = (0.5/2700)*t_m1 + 79290/2700;
else
    T6 = 30.7;
end

% Layer 2
if t_m1>=0 && t_m1<1800
    L2 = (11.2/1800)*t_m1 + (36720/1800);
elseif t_m1>=1800 && t_m1<4500
    L2 = (1.3/2700)*t_m1 + 82980/2700;
elseif t_m1>=4500 && t_m1<7200
    L2 = (1.7/2700)*t_m1 + 81180/2700;
else
    L2 = 34.6;
end

% Layer 1
if t_m1>=0 && t_m1<1800
    L1 = (13.4/1800)*t_m1 + (36720/1800);
elseif t_m1>=1800 && t_m1<4500
    L1 = (1/2700)*t_m1 + 89460/2700;
elseif t_m1>=4500 && t_m1<7200
    L1 = (1.3/2700)*t_m1 + 88110/2700;
else
    L1 = 36.1;
end

T = zeros(14,9,3);
T(:,:,1) = L1;
T(:,:,2) = L2;
T(1:4,1:5,3) = T1; T(1:4,6:9,3) = T2;
T(5:9,1:5,3) = T3; T(5:9,6:9,3) = T4;
T(10:14,1:5,3) = T5; T(10:14,6:9,3) = T6;

T = T + 273; % Temperature conversion from celsius to kelvin 
T_in = T(1,3,2); % Approximate inlet gas temperature to chamber (sourced from gas tanks)
%% Define initial condition for temperature-dependent parameters
Temp = [25+273 50+273]; 
% Diffusion coefficient of oxygen in gas corresponding to Temp
D_range = [0.205 0.237]; % Calculated using Hirschfelder equation
D_int = griddedInterpolant(Temp,D_range,'linear');

for z=1:3
    for y=1:14
        for x=1:9
            D(y,x,z) = D_int(T(y,x,z));
        end
    end
end

%% Flow Profile Calculation
% Load file with calculated flow splits from ANSYS data
load("C:\Users\Nida Qayyum\Documents\MATLAB\Flow_data.mat",'x_y_por_l2','x_y_por_l3','x_n_avg_l3_2','x_p_avg_l3_2','x_n_avg','x_p_avg','x_p_avg_l1','x_n_avg_l1','x_p_avg_l3','x_n_avg_l3')

Qin = 450/60; % Total inlet flow to chamber
Q = zeros(14,9,3);
Inlet_l1 = 1-Inlet_l2-Inlet_l3; % Fraction of flow to Layer 1

% Flow profile for Layer 2
z = 2;
for y=1:14
    if y == 1
        Q(y,3,z) = Qin; % Total inlet to chamber
        Q(y,2,z) = Qin*Inlet_l2*(abs(x_n_avg)/(abs(x_n_avg) + abs(x_p_avg)));
        Q(y,4,z) = Qin*Inlet_l2*(abs(x_p_avg)/(abs(x_n_avg) + abs(x_p_avg)));
        Q(y,1,z) = Q(y,2,z);
        Q(y,5,z) = Q(y,4,z);
        Q(y,6,z) = Q(y,5,z)*x_y_por_l2(y,5);
        Q(y,7,z) = Q(y,6,z);
        Q(y,8,z) = Q(y,7,z);
        Q(y,9,z) = Q(y,8,z);
    elseif y == 2
        for x = 1:9
            if x == 1 || x == 9
                Q(y,x,z) = Q(y-1,x,z);
            elseif x == 5
                Q(y,x,z) = Q(y-1,x,z)*(1 - x_y_por_l2(y-1,x));
            else
                Q(y,x,z) = 0;
            end
        end
    elseif y == 3 || y == 4 || y == 7 || y == 8 || y == 12 || y == 13
        for x = 1:9
            if x == 1 || x == 5 || x == 9
                Q(y,x,z) = Q(y-1,x,z);
            else
                Q(y,x,z) = 0;
            end
        end
    elseif y == 5
        Q(y,1,z) = Q(y-1,1,z);
        Q(y,2,z) = Q(y,1,z)*x_y_por_l2(y,1);
        Q(y,3,z) = Q(y,2,z);
        Q(y,4,z) = Q(y,3,z);
        Q(y,5,z) = Q(y-1,5,z) + Q(y,4,z);
        Q(y,6,z) = Q(y,5,z)*x_y_por_l2(y,5);
        Q(y,7,z) = Q(y,6,z);
        Q(y,8,z) = Q(y,7,z);
        Q(y,9,z) = Q(y-1,9,z) + Q(y,8,z);
    elseif y == 6 || y == 11
        for x = 1:9
            if x == 9
                Q(y,x,z) = Q(y-1,x,z);
            elseif x == 1 || x == 5
                Q(y,x,z) = Q(y-1,x,z)*(1-x_y_por_l2(y-1,x));
            else
                Q(y,x,z) = 0;
            end
        end
    elseif y == 9
        Q(y,1,z) = Q(y-1,1,z);
        Q(y,2,z) = Q(y,1,z)*x_y_por_l2(y,1);
        Q(y,3,z) = Q(y,2,z)*x_y_por_l2(y,2);
        Q(y,4,z) = Q(y,3,z)*x_y_por_l2(y,3);
        Q(y,5,z) = Q(y-1,5,z) + Q(y,4,z)*x_y_por_l2(y,4);
        Q(y,6,z) = Q(y,5,z)*x_y_por_l2(y,5);
        Q(y,7,z) = Q(y,6,z)*x_y_por_l2(y,6);
        Q(y,8,z) = Q(y,7,z)*x_y_por_l2(y,7);
        Q(y,9,z) = Q(y-1,9,z) + Q(y,8,z)*x_y_por_l2(y,8);
    elseif y == 10
        for x = 1:9
            if x == 1
                Q(y,x,z) = Q(y-1,x,z)*(1 - x_y_por_l2(y-1,x));
            elseif x == 2 || x == 6
                Q(y,x,z) = Q(y,x-1,z)*x_y_por_l2(y,x-1) + Q(y-1,x,z)*(1 - x_y_por_l2(y-1,x));
            elseif x == 9
                Q(y,x,z) = Q(y-1,x,z) + Q(y,x-1,z);
            else
                Q(y,x,z) = Q(y,x-1,z) + Q(y-1,x,z)*(1 - x_y_por_l2(y-1,x));
            end
        end
    else
        Q(y,1,z) = Q(y-1,1,z);
        Q(y,2,z) = Q(y,1,z);
        Q(y,3,z) = Q(y,2,z);
        Q(y,4,z) = Q(y,3,z);
        Q(y,5,z) = Q(y,4,z) + Q(y-1,5,z);
        Q(y,6,z) = Q(y,5,z);
        Q(y,9,z) = Q(y-1,9,z);
        Q(y,8,z) = Q(y,9,z);
        Q(y,7,z) = Q(y,6,z) + Q(y,8,z); % Adjust with flow rates from top and bottom after defining Layer 1 and 3 profile
    end
end

% Flow profile for Layer 1
z = 1;
for y=1:14
    if y == 1
        Q(y,3,z) = Qin*Inlet_l1; % Total inlet to Layer 1
        Q(y,2,z) = Qin*Inlet_l1*(abs(x_n_avg_l1)/(abs(x_n_avg_l1) + abs(x_p_avg_l1)));
        Q(y,4,z) = Qin*Inlet_l1*(abs(x_p_avg_l1)/(abs(x_n_avg_l1) + abs(x_p_avg_l1)));
        Q(y,1,z) = Q(y,2,z);
        Q(y,5,z) = Q(y,4,z);
        Q(y,6,z) = Q(y,5,z);
        Q(y,7,z) = Q(y,6,z);
        Q(y,8,z) = Q(y,7,z);
        Q(y,9,z) = Q(y,8,z);
    elseif y>1 && y<14
        for x = 1:9
            if x == 1 || x == 9
                Q(y,x,z) = Q(y-1,x,z);
            else
                Q(y,x,z) = 0;
            end
        end
    else
        Q(y,1,z) = Q(y-1,1,z);
        Q(y,2,z) = Q(y,1,z);
        Q(y,3,z) = Q(y,2,z);
        Q(y,4,z) = Q(y,3,z);
        Q(y,5,z) = Q(y,4,z);
        Q(y,6,z) = Q(y,5,z);
        Q(y,9,z) = Q(y-1,9,z);
        Q(y,8,z) = Q(y,9,z);
        Q(y,7,z) = Q(y,6,z) + Q(y,8,z); 
    end
end

% Flow profile for Layer 3
z = 3;
for y=1:14
    if y == 1
        Q(y,3,z) = Qin*Inlet_l3; % Total inlet to Layer 3
        Q(y,2,z) = Qin*Inlet_l3*x_y_por_l3(y,3)*(abs(x_n_avg_l3)/(abs(x_n_avg_l3) + abs(x_p_avg_l3)));
        Q(y,4,z) = Qin*Inlet_l3*x_y_por_l3(y,3)*(abs(x_p_avg_l3)/(abs(x_n_avg_l3) + abs(x_p_avg_l3)));
        Q(y,1,z) = Q(y,2,z)*x_y_por_l3(y,2);
        Q(y,5,z) = Q(y,4,z)*x_y_por_l3(y,4);
        Q(y,6,z) = Q(y,5,z)*x_y_por_l3(y,5);
        Q(y,7,z) = Q(y,6,z)*x_y_por_l3(y,6);
        Q(y,8,z) = Q(y,7,z)*x_y_por_l3(y,7);
        Q(y,9,z) = Q(y,8,z)*x_y_por_l3(y,8);
    elseif y == 2
        Q(y,3,z) = Q(y-1,3,z)*(1 - x_y_por_l3(y-1,3));
        Q(y,2,z) = Q(y,3,z)*x_y_por_l3(y,3)*(abs(x_n_avg_l3_2)/(abs(x_n_avg_l3_2) + abs(x_p_avg_l3_2))) + Q(y-1,2,z)*(1 - x_y_por_l3(y-1,2));
        Q(y,4,z) = Q(y,3,z)*x_y_por_l3(y,3)*(abs(x_p_avg_l3_2)/(abs(x_n_avg_l3_2) + abs(x_p_avg_l3_2))) + Q(y-1,4,z)*(1 - x_y_por_l3(y-1,4));
        Q(y,1,z) = Q(y-1,1,z) + Q(y,2,z)*x_y_por_l3(y,2);
        Q(y,5,z) = Q(y,4,z)*x_y_por_l3(y,4) + Q(y-1,5,z)*(1 - x_y_por_l3(y-1,5));
        Q(y,6,z) = Q(y,5,z)*x_y_por_l3(y,5) + Q(y-1,6,z)*(1 - x_y_por_l3(y-1,6));
        Q(y,7,z) = Q(y,6,z)*x_y_por_l3(y,6) + Q(y-1,7,z)*(1 - x_y_por_l3(y-1,7));
        Q(y,8,z) = Q(y,7,z)*x_y_por_l3(y,7) + Q(y-1,8,z)*(1 - x_y_por_l3(y-1,8));
        Q(y,9,z) = Q(y-1,9,z) +  Q(y,8,z)*x_y_por_l3(y,8);
    elseif y == 3
        Q(y,1,z) = Q(y-1,1,z);
        Q(y,2,z) = Q(y,1,z)*x_y_por_l3(y,1) + Q(y-1,2,z)*(1 - x_y_por_l3(y-1,2));
        Q(y,3,z) = Q(y,2,z)*x_y_por_l3(y,2) + Q(y-1,3,z)*(1 - x_y_por_l3(y-1,3));
        Q(y,4,z) = Q(y,3,z)*x_y_por_l3(y,3) + Q(y-1,4,z)*(1 - x_y_por_l3(y-1,4));
        Q(y,5,z) = Q(y,4,z)*x_y_por_l3(y,4) + Q(y-1,5,z)*(1 - x_y_por_l3(y-1,5));
        Q(y,6,z) = Q(y,5,z)*x_y_por_l3(y,5) + Q(y-1,6,z)*(1 - x_y_por_l3(y-1,6));
        Q(y,9,z) = Q(y-1,9,z);
        Q(y,8,z) = Q(y,9,z)*x_y_por_l3(y,9) + Q(y-1,8,z)*(1 - x_y_por_l3(y-1,8));
        Q(y,7,z) = Q(y,6,z)*x_y_por_l3(y,6) + Q(y-1,7,z)*(1 - x_y_por_l3(y-1,7)) + Q(y,8,z)*x_y_por_l3(y,8);
    elseif y>3 && y<14
        Q(y,1,z) = Q(y-1,1,z)*(1 - x_y_por_l3(y-1,1));
        for x=2:6
            Q(y,x,z) = Q(y,x-1,z)*x_y_por_l3(y,x-1) + Q(y-1,x,z)*(1 - x_y_por_l3(y-1,x));
        end
        Q(y,9,z) = Q(y-1,9,z)*(1 - x_y_por_l3(y-1,9));
        Q(y,8,z) = Q(y,9,z)*x_y_por_l3(y,9) + Q(y-1,8,z)*(1 - x_y_por_l3(y-1,8));
        Q(y,7,z) = Q(y-1,7,z) + Q(y,8,z)*x_y_por_l3(y,8) + Q(y,6,z)*x_y_por_l3(y,6);
    else
        Q(y,1,z) = Q(y-1,1,z)*(1 - x_y_por_l3(y-1,1));
        Q(y,2,z) = Q(y,1,z) + Q(y-1,2,z)*(1 - x_y_por_l3(y-1,2));
        Q(y,3,z) = Q(y,2,z) + Q(y-1,3,z)*(1 - x_y_por_l3(y-1,3));
        Q(y,4,z) = Q(y,3,z) + Q(y-1,4,z)*(1 - x_y_por_l3(y-1,4));
        Q(y,5,z) = Q(y,4,z) + Q(y-1,5,z)*(1 - x_y_por_l3(y-1,5));
        Q(y,6,z) = Q(y,5,z) + Q(y-1,6,z)*(1 - x_y_por_l3(y-1,6));
        Q(y,9,z) = Q(y-1,9,z)*(1 - x_y_por_l3(y-1,9));
        Q(y,8,z) = Q(y,9,z) + Q(y-1,8,z)*(1 - x_y_por_l3(y-1,8));
        Q(y,7,z) = Q(y-1,7,z) + Q(y,6,z) + Q(y,8,z);
    end
end

Q(14,7,2) = Q(14,7,2) + Q(14,7,1) + Q(14,7,3); % Total outlet from chamber

%% Define initial conditions, non temperature-dependent parameters, and solving parameters
% Heights for Layers 1-3
h = zeros([14,9,3]);
h(:,:,1) = 1;
h(:,:,2) = 1.3; 
h(:,:,3) = 0.4; 

% Parameters
P = 760; % Total air pressure (mmHg)
R = 62360; % Ideal gas constant (mL mmHg/mol K)

% Solving parameters
%tend = 120; % End-time for step change 
tend = tf(end); % End-time for IH
dt = 0.01; % Model time step
s = round(tend/dt);
dx = 1; dy = 1;
t = zeros(1,s);
%protocol = 0; % For step change
protocol = 1; % For IH

c = zeros([14,9,3,s]); % Oxygen concentration in gas phase
f0_in = 0.21; % Inlet oxygen fraction
cO2in(1) = f0_in*P/(R*T_in);


%% Solving loop for time-dependent oxygen concentration profile in the chamber
for i=2:s
    t(i) = t(i-1) + dt;
    
    %% Temperature Data
    t_m1 = t_m1 + dt; t_curve = t_m1; % profile time to be used with measured temperature data

    % Well 1
    if t_curve>=0 && t_curve<1800
        T(1:4,1:5,3) = (6.9/1800)*t_curve + 20.3;
    elseif t_curve>=1800 && t_curve<4500
        T(1:4,1:5,3) = (1.6/2700)*t_curve + (70560/2700);
    elseif t_curve>=4500 && t_curve<7200
        T(1:4,1:5,3) = (0.5/2700)*t_curve + (75510/2700);
    else
        T(1:4,1:5,3) = 29.3;
    end

    % Well 2
    if t_curve>=0 && t_curve<1800
        T(1:4,6:9,3) = (7.6/1800)*t_curve + 36720/1800;
    elseif t_curve>=1800 && t_curve<4500
        T(1:4,6:9,3) = (1.3/2700)*t_curve + 73260/2700;
    elseif t_curve>=4500 && t_curve<7200
        T(1:4,6:9,3) = (0.1/2700)*t_curve + 78660/2700;
    else
        T(1:4,6:9,3) = 29.4;
    end


    % Well 3
    if t_curve>=0 && t_curve<1800
        T(5:9,1:5,3) = (8.6/1800)*t_curve + 36900/1800;
    elseif t_curve>=1800 && t_curve<4500
        T(5:9,1:5,3) = (1.6/2700)*t_curve + 75690/2700;
    elseif t_curve>=4500 && t_curve<7200
        T(5:9,1:5,3) = -(0.1/2700)*t_curve + 83340/2700;
    else
        T(5:9,1:5,3) = 30.6;
    end

    % Well 4
    if t_curve>=0 && t_curve<1800
        T(5:9,6:9,3) = (9.3/1800)*t_curve + 36720/1800;
    elseif t_curve>=1800 && t_curve<4500
        T(5:9,6:9,3) = (1.1/2700)*t_curve + 78210/2700;
    elseif t_curve>=4500 && t_curve<7200
        T(5:9,6:9,3) = (0.4/2700)*t_curve + 81360/2700;
    else
        T(5:9,6:9,3) = 31.2;
    end

    % Well 5
    if t_curve>=0 && t_curve<1800
        T(10:14,1:5,3) = (8.5/1800)*t_curve + 20.4;
    elseif t_curve>=1800 && t_curve<4500
        T(10:14,1:5,3) = (1.1/2700)*t_curve + 76050/2700;
    elseif t_curve>=4500 && t_curve<7200
        T(10:14,1:5,3) = (0.5/2700)*t_curve + 78750/2700;
    else
        T(10:14,1:5,3) = 30.5;
    end

    % Well 6
    if t_curve>=0 && t_curve<1800
        T(10:14,6:9,3) = (8.6/1800)*t_curve + 20.5;
    elseif t_curve>=1800 && t_curve<4500
        T(10:14,6:9,3) = (1.1/2700)*t_curve + 76590/2700;
    elseif t_curve>=4500 && t_curve<7200
        T(10:14,6:9,3) = (0.5/2700)*t_curve + 79290/2700;
    else
        T(10:14,6:9,3) = 30.7;
    end

    % Layer 2
    if t_curve>=0 && t_curve<1800
        T(:,:,2) = (11.2/1800)*t_curve + (36720/1800);
    elseif t_curve>=1800 && t_curve<4500
        T(:,:,2) = (1.3/2700)*t_curve + 82980/2700;
    elseif t_curve>=4500 && t_curve<7200
        T(:,:,2) = (1.7/2700)*t_curve + 81180/2700;
    else
        T(:,:,2) = 34.6;
    end
    
    % Layer 1
    if t_curve>=0 && t_curve<1800
        T(:,:,1) = (13.4/1800)*t_curve + (36720/1800);
    elseif t_curve>=1800 && t_curve<4500
        T(:,:,1) = (1/2700)*t_curve + 89460/2700;
    elseif t_curve>=4500 && t_curve<7200
        T(:,:,1) = (1.3/2700)*t_curve + 88110/2700;
    else
        T(:,:,1) = 36.1;
    end

    T = T + 273;

    %% Temperature-dependent parameter values
    D1 = D_int(T(2,2,3)); D(1:4,1:5,3) = D1; 
    D2 = D_int(T(2,6,3)); D(1:4,6:9,3) = D2;
    D3 = D_int(T(6,2,3)); D(5:9,1:5,3) = D3; 
    D4 = D_int(T(6,6,3)); D(5:9,6:9,3) = D4; 
    D5 = D_int(T(11,2,3)); D(10:14,1:5,3) = D5;
    D6 = D_int(T(11,6,3)); D(10:14,6:9,3) = D6;
    DL1 = D_int(T(1,1,1)); D(:,:,1) = DL1;
    DL2 = D_int(T(1,1,2)); D(:,:,2) = DL2;

    %% Layer 2
    z = 2;
    for y=1:14
        if y == 1
            for x=1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dx*dy)*(c(y,x,z-1,i-1)-c(y,x,z,i-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + D(y,x,z)*(dx*dy)*(c(y,x,z+1,i-1)-c(y,x,z,i-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + c(y,x,z,i-1);
                elseif x == 2
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,i-1)-c(y,x,z,i-1))/dx + D(y,x,z)*(dx*dy)*(c(y,x,z-1,i-1)-c(y,x,z,i-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + D(y,x,z)*(dx*dy)*(c(y,x,z+1,i-1)-c(y,x,z,i-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + c(y,x,z,i-1);
                elseif x == 3
                    if protocol == 1
                        c(y,x,z,i) = interp1(tf,O2_f*(P/(R*T_in)),t(i),'next');
                    else
                        c(y,x,z,i) = cO2in;
                    end
                elseif x>3 && x<5 || x>5 && x<9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx + D(y,x,z)*(dx*dy)*(c(y,x,z-1,i-1)-c(y,x,z,i-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + D(y,x,z)*(dx*dy)*(c(y,x,z+1,i-1)-c(y,x,z,i-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + c(y,x,z,i-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dx*dy)*(c(y,x,z-1,i-1)-c(y,x,z,i-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + D(y,x,z)*(dx*dy)*(c(y,x,z+1,i-1)-c(y,x,z,i-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + c(y,x,z,i-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dx*dy)*(c(y,x,z-1,i-1)-c(y,x,z,i-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + D(y,x,z)*(dx*dy)*(c(y,x,z+1,i-1)-c(y,x,z,i-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + c(y,x,z,i-1);
                end
            end
        elseif y == 2
            for x=1:9
                if x == 1 || x == 9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dx*dy)*(c(y,x,z-1,i-1)-c(y,x,z,i-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + D(y,x,z)*(dx*dy)*(c(y,x,z+1,i-1)-c(y,x,z,i-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + c(y,x,z,i-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dx*dy)*(c(y,x,z+1,i-1)-c(y,x,z,i-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + c(y,x,z,i-1);
                else
                    c(y,x,z,i) = 0;
                end
            end
        elseif y>2 && y<5 || y>5 && y<9 || y>10 && y<14
            for x=1:9
                if x == 1 || x == 5 || x == 9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*dt/V + c(y,x,z,i-1);
                else
                    c(y,x,z,i) = 0;
                end
            end
        elseif y == 5 
            for x = 1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx + D(y,x,z)*(dx*(h(y,x,z)))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1);
                elseif x>1 && x<5 || x>5 && x<9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x-1,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx + Q(y-1,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy + Q(y,x-1,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1);
                end
            end
        elseif y == 9 
            for x = 1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1);
                elseif x>1 && x<5 || x>5 && x<9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx + D(y,x,z)*(dy*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x-1,z)*x_y_por_l2(y,x-1)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx + Q(y-1,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1); 
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy + Q(y,x-1,z)*x_y_por_l2(y,x-1)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1);
                end
            end
        elseif y == 10
            for x = 1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1);
                elseif x == 2 || x == 6
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x-1,z)*x_y_por_l2(y,x-1)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx + Q(y-1,x,z)*(1 - x_y_por_l2(y-1,x))*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x-1,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx + Q(y-1,x,z)*(1 - x_y_por_l2(y-1,x))*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1);
                elseif x == 9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy + Q(y,x-1,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x-1,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx + Q(y-1,x,z)*(1 - x_y_por_l2(y-1,x))*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                end
            end
        else
            for x = 1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                elseif x>1 && x<5 || x == 6
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x-1,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx + Q(y-1,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x-1,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx + Q(y,x+1,z)*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx + Q(y,x,z-1)*c(y,x,z-1,i-1) - D(y,x,z)*(dx*dy)*(c(y,x,z,i-1)-c(y,x,z-1,i-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + Q(y,x,z+1)*c(y,x,z+1,i-1) - D(y,x,z)*(dx*dy)*(c(y,x,z,i-1)-c(y,x,z+1,i-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*c(y,x,z,i-1))*(dt/V) + c(y,x,z,i-1); %one outlet, multiple inlets
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                end
            end
        end
    end

    %% Layer 1
    z = 1;

    for y=1:14
        if y == 1
            for x=1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx - D(y,x,z)*(dx*dy)*(c(y,x,z,i-1)-c(y,x,z+1,i-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1);
                elseif x == 2
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx - D(y,x,z)*(dx*dy)*(c(y,x,z,i-1)-c(y,x,z+1,i-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                elseif x == 3
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x,z+1,i-1) - D(y,x,z)*(dx*dy)*(c(y,x,z,i-1)-c(y,x,z+1,i-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,i-1)-c(y,x,z,i-1))/dx + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1); %inlet from top
                elseif x>3 && x<9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - D(y,x,z)*(dx*dy)*(c(y,x,z,i-1)-c(y,x,z+1,i-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - D(y,x,z)*(dx*dy)*(c(y,x,z,i-1)-c(y,x,z+1,i-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1);
                end
            end
        elseif y == 2
            for x=1:9
                if x == 1 || x == 9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - D(y,x,z)*(dx*dy)*(c(y,x,z,i-1)-c(y,x,z+1,i-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1);
                else
                    c(y,x,z,i) = 0;
                end
            end
        elseif y>2 && y<14
            for x=1:9
                if x == 1 || x == 9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1);
                else
                    c(y,x,z,i) = 0;
                end
            end
        else
            for x = 1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                elseif x>1 && x<7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x-1,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx + Q(y,x+1,z)*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*dy)*(c(y,x,z+1,i-1)-c(y,x,z,i-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*(dt/V) + c(y,x,z,i-1); %one outlet, multiple inlets
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                end
            end
        end
    end


    %% Layer 3
    z = 3;
 
    for y=1:14
        if y == 1
            for x=1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx - D(y,x,z)*(dx*dy)*(c(y,x,z,i-1)-c(y,x,z-1,i-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1); 
                elseif x == 2
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx - D(y,x,z)*(dx*dy)*(c(y,x,z,i-1)-c(y,x,z-1,i-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                elseif x == 3
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x,z-1,i-1) - D(y,x,z)*(dx*dy)*(c(y,x,z,i-1)-c(y,x,z-1,i-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,i-1)-c(y,x,z,i-1))/dx + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1); 
                elseif x>3 && x<9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - D(y,x,z)*(dx*dy)*(c(y,x,z,i-1)-c(y,x,z-1,i-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - D(y,x,z)*(dx*dy)*(c(y,x,z,i-1)-c(y,x,z-1,i-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1); 
                end
            end
        elseif y == 2
            for x=1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - D(y,x,z)*(dx*dy)*(c(y,x,z,i-1)-c(y,x,z-1,i-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + Q(y,x+1,z)*x_y_por_l3(y,x+1)*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1); 
                elseif x == 2
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*(abs(x_n_avg_l3_2)/(abs(x_n_avg_l3_2) + abs(x_p_avg_l3_2)))*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                elseif x == 3
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,i-1)-c(y,x,z,i-1))/dx + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                elseif x == 4
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*(abs(x_p_avg_l3_2)/(abs(x_n_avg_l3_2) + abs(x_p_avg_l3_2)))*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - D(y,x,z)*(dx*dy)*(c(y,x,z,i-1)-c(y,x,z-1,i-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + Q(y,x-1,z)*x_y_por_l3(y,x-1)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                elseif x>5 && x<9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - D(y,x,z)*(dx*dy)*(c(y,x,z,i-1)-c(y,x,z-1,i-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + Q(y,x-1,z)*x_y_por_l3(y,x-1)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1);
                end
            end
        elseif y == 3
            for x=1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1); 
                elseif x>1 && x<7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx + Q(y,x+1,z)*x_y_por_l3(y,x+1)*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1);
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1); 
                end
            end
        elseif y>3 && y<14
            for x=1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1); 
                elseif x>1 && x<7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx + Q(y,x+1,z)*x_y_por_l3(y,x+1)*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy)*(dt/V) + c(y,x,z,i-1);
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,i-1)-c(y,x,z,i-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1); 
                end
            end
        else
            for x=1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1); 
                elseif x>1 && x<7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy + Q(y,x-1,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy + Q(y,x-1,z)*c(y,x-1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x-1,z,i-1))/dx + Q(y,x+1,z)*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dx*dy)*(c(y,x,z-1,i-1)-c(y,x,z,i-1))/(0.5*(h(y,x,z)+h(y,x,z-1))))*(dt/V) + c(y,x,z,i-1);
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy + Q(y,x+1,z)*c(y,x+1,z,i-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,i-1)-c(y,x+1,z,i-1))/dx - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,i) = (Q(y,x,z)*c(y-1,x,z,i-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,i-1)-c(y-1,x,z,i-1))/dy - Q(y,x,z)*c(y,x,z,i-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,i-1)-c(y,x,z,i-1))/dx)*(dt/V) + c(y,x,z,i-1); 
                end
            end
        end
    end
end

%% Outputs
TR_p = 0.5*(squeeze(c(2,8,3,:))*R*T(2,8,3)*100/P + squeeze(c(2,9,3,:))*R*T(2,9,3)*100/P);
BL_p = 0.5*(squeeze(c(13,1,3,:))*R*T(13,1,3)*100/P + squeeze(c(13,2,3,:))*R*T(13,2,3)*100/P);
BR_p = 0.5*(squeeze(c(13,6,3,:))*R*T(13,6,3)*100/P + squeeze(c(13,7,3,:))*R*T(13,7,3)*100/P);
TR_10 = TR_p(round(10/dt)-round(td/dt)); TR_30 = TR_p(round(30/dt)-round(td/dt)); TR_60 = TR_p(round(60/dt)-round(td/dt)); TR_90 = TR_p(round(90/dt)-round(td/dt));
BL_10 = BL_p(round(10/dt)-round(td/dt)); 
BL_30 = BL_p(round(30/dt)-round(td/dt)); BL_60 = BL_p(round(60/dt)-round(td/dt)); BL_90 = BL_p(round(90/dt)-round(td/dt));
%hold on
%plot(t+td-371.3,BL_p,'Color',cp)
%end