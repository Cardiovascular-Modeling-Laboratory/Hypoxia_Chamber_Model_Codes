%% Prepare parameter set for model run (FIRST RUN CHAMBER_LIQUID_PARAMETERS_FILTER_PAPER.M)
% Using sample parameter set from filtered list
for i=1:1
    row = pars_run(i);
    T1 = x_pars_v(row,1);
    k1 = x_pars_v(row,2);
    k2 = x_pars_v(row,3);
    k3 = x_pars_v(row,4);
    k4 = x_pars_v(row,5);
    k5 = x_pars_v(row,6);
    k6 = x_pars_v(row,7);
    Inlet_l2 = x_pars_v(row,8);
    Inlet_l3 = x_pars_v(row,9);
end

%% Temperature data
% Find time on temperature curve using Well 1
if T1>=20.7 && T1<23.4
    t_m1 = (2000/3)*T1 - 13800;
elseif T1>=23.4 && T1<25.1
    t_m1 = (27000/17)*T1 - 601200/17;
elseif T1>=25.1 && T1<25.5
    t_m1 = 6750*T1 - 164925;
else
    t_m1 = 7200;
end

% Determine gas temperatures above remaining wells and for Layers 1 and 2
% Well 2
if t_m1>=0 && t_m1<1800
    T2 = (7/3000)*t_m1 + 104/5;
elseif t_m1>=1800 && t_m1<4500
    T2 = (7/27000)*t_m1 + 368/15;
elseif t_m1>=4500 && t_m1<7200
    T2 = (1/5400)*t_m1 + 373/15;
else
    T2 = 26.2;
end

% Well 3
if t_m1>=0 && t_m1<1800
    T3 = (83/18000)*t_m1 + 103/5;
elseif t_m1>=1800 && t_m1<4500
    T3 = (17/27000)*t_m1 + 833/30;
elseif t_m1>=4500 && t_m1<7200
    T3 = -(0.2/2700)*t_m1 + 83520/2700;
else
    T3 = 30.4;
end

% Well 4
if t_m1>=0 && t_m1<1800
    T4 = (31/6000)*t_m1 + 103/5;
elseif t_m1>=1800 && t_m1<4500
    T4 = (7/13500)*t_m1 + 869/30;
else
    T4 = 31.3;
end

% Well 5
if t_m1>=0 && t_m1<1800
    T5 = (41/9000)*t_m1 + 20.7;
elseif t_m1>=1800 && t_m1<4500
    T5 = (1/1350)*t_m1 + 827/30;
elseif t_m1>=4500 && t_m1<7200
    T5 = (1/6750)*t_m1 + 907/30;
else
    T5 = 31.3;
end

% Well 6
if t_m1>=0 && t_m1<1800
    T6 = (7/1500)*t_m1 + 21;
elseif t_m1>=1800 && t_m1<4500
    T6 = (17/27000)*t_m1 + 424/15;
elseif t_m1>=4500 && t_m1<7200
    T6 = (1/27000)*t_m1 + 464/15;
else
    T6 = 31.2;
end

% Layer 2
if t_m1>=0 && t_m1<1800
    L2 = (7.6/1800)*t_m1 + 37800/1800;
elseif t_m1>=1800 && t_m1<4500
    L2 = (2.5/2700)*t_m1 + 72720/2700;
elseif t_m1>=4500 && t_m1<7200
    L2 = (0.5/2700)*t_m1 + 81720/2700;
else
    L2 = 31.6;
end

% Layer 1
if t_m1>=0 && t_m1<1800
    L1 = (10.4/1800)*t_m1 + 37440/1800;
elseif t_m1>=1800 && t_m1<4500
    L1 = (1.2/2700)*t_m1 + 82080/2700;
elseif t_m1>=4500 && t_m1<7200
    L1 = (1.1/2700)*t_m1 + 82530/2700;
else
    L1 = 33.5;
end

% Gas temperatures
T = zeros(14,9,3);
T(:,:,1) = L1;
T(:,:,2) = L2;
T(1:4,1:5,3) = T1; T(1:4,6:9,3) = T2;
T(5:9,1:5,3) = T3; T(5:9,6:9,3) = T4;
T(10:14,1:5,3) = T5; T(10:14,6:9,3) = T6;
T = T + 273; % Temperature conversion from celsius to kelvin
T_in = T(1,3,2); % Approximate inlet gas temperature to chamber (sourced from gas tanks)

% Liquid temperatures inside wells (determined from regression of measured gas vs. liquid temeperatures)
T_w = zeros(14,9,3);
T_w(:,:,1) = 0; % No liquid in Layer 1
T_w(:,:,2) = 0; % Initialize
T_w(2:4,2:4,2) = 1.1*T1 - 2.11; T_w(2:4,6:8,2) = 1.1*T2 - 2.11;
T_w(6:8,2:4,2) = 1.1*T3 - 2.11; T_w(6:8,6:8,2) = 1.1*T4 - 2.11;
T_w(11:13,2:4,2) = 1.1*T5 - 2.11; T_w(11:13,6:8,2) = 1.1*T6 - 2.11;
T_w(:,:,3) = 0; % No liquid in Layer 3
T_w = T_w + 273; %temperature conversion from celsius to kelvin

% Chamber lid temperatures (to be determined using Layer 3 temperatures and experimental data of difference factor)
T_lid = T(:,:,3); T_diff = zeros(14,9,1); % Initialize

%% Define initial condition for temperature-dependent parameters
Temp = [25+273 50+273];
% Diffusion coefficient for oxygen in gas corresponding to Temp
D_range = [0.205 0.237]; % Calculated using Hirschfelder equation
% Saturation pressure of water vapor corresponding to Temp
Psat_range = [23.8 92.6];
% Diffusion coefficient for water vapor in gas corresponding to Temp
Dw_range = [0.26 0.303]; % Calculated using Hirschfelder equation

% Diffusion coefficient of oxygen in water
Temp_DO2 = [25+273 40+273];
% Diffusion coefficient for oxygen in water corresponding to Temp_DO2
DO2_range = [2.2e-5 3.3e-5]; % (Source: Oxygen Diffusion in Water. Estimation of the Coefficients  of Molecular Diffusion -- AIP Publishing)

% Solubility of oxygen in water
Temp_sol = [25+273 35+273];
% Solubulity of oxygen in water corresponding to Temp_sol
sol_range = [259 218]; %umol/kg (Source: The solubility of oxygen in the major sea salts and their mixtures at 25°C)
Beta_O2_range = [sol_range(1)*1e-9/(0.21*(760-interp1(Temp,Psat_range,Temp_sol(1),'linear','extrap'))) sol_range(2)*1e-9/(0.21*(760-interp1(Temp,Psat_range,Temp_sol(2),'linear','extrap')))]; %units are mol/mL mmHg

D = zeros(14,9,3); Psat = zeros(14,9,3); Dw = zeros(14,9,3); DO2 = zeros(14,9,3); Beta_O2 = zeros(14,9,3); Psat_cond = zeros(14,9,3);
D_int = griddedInterpolant(Temp,D_range,'linear');
Psat_int = griddedInterpolant(Temp,Psat_range,'linear');
Dw_int = griddedInterpolant(Temp,Dw_range,'linear');
DO2_int = griddedInterpolant(Temp_DO2,DO2_range,'linear');
Beta_int = griddedInterpolant(Temp_sol,Beta_O2_range,'linear');

for z=1:3
    for y=1:14
        for x=1:9
            D(y,x,z) = interp1(Temp,D_range,T(y,x,z),'linear','extrap'); % Diffusion of oxygen in gas phase
            if T_w(y,x,z) == 273 % No liquid present
                Psat(y,x,z) = interp1(Temp,Psat_range,T(y,x,z),'linear','extrap'); % Use gas temperature for saturation pressure
            else
                Psat(y,x,z) = interp1(Temp,Psat_range,T_w(y,x,z),'linear','extrap'); % Use liquid temperature for saturation pressure of water at surface
            end
            Dw(y,x,z) = interp1(Temp,Dw_range,T(y,x,z),'linear','extrap'); % Diffusion of water vapor in gas phase
            DO2(y,x,z) = interp1(Temp_DO2,DO2_range,T_w(y,x,z),'linear','extrap'); % Diffusion of oxygen in water
            Beta_O2(y,x,z) = interp1(Temp_sol,Beta_O2_range,T_w(y,x,z),'linear','extrap'); % Solubility of oxygen in water
        end
    end
end
D1_O2 = DO2(2,2,2); D2_O2 = DO2(2,6,2); D3_O2 = DO2(6,2,2); D4_O2 = DO2(6,6,2); D5_O2 = DO2(11,2,2); D6_O2 = DO2(11,6,2);

%% Flow Profile Calculation
% Load file with calculated flow splits from ANSYS data
[flow_file, flow_path] = uigetfile('*.mat','Select file with results of velocity_data_extraction_PAPER (titled Flow_data.mat)');
load([flow_path flow_file],'x_y_por_l2','x_y_por_l3','x_n_avg_l3_2','x_p_avg_l3_2','x_n_avg','x_p_avg','x_p_avg_l1','x_n_avg_l1','x_p_avg_l3','x_n_avg_l3')

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

% Well dimensions
P_h = 0.9; % Height of petri dishes (wells)
S_w = 9; % Model surface area approximation of wells

% Parameters
P = 760; % Total air pressure (mmHg)
R = 62360; % Ideal gas constant (mL mmHg/mol K)
Mw = 18; % Molecular weight of water (g/mol)
rho = 1; % Density of water (g/mL)

% Cellular metabolics (using Michaelis-Menten method)
No_Cells = 650000; % Number of cells per well
O2Met = 8e-17; % Oxygen consumption per cell (mol/cell s)
O2consump = O2Met*No_Cells; % Total oxygen consumption per well
ke = 5.5e-10; % Michaelis-Menten constant (mol/cm3)


% Solving parameters
protocol = input('Are you running a step change (=0) or IH protocol (=1)?');
if protocol == 0
    tend = input('Enter total time of run in seconds:');
    f0_in = input('Enter the inlet oxygen fraction:');
else
    [filename_IH, path_IH] = uigetfile('*.*','Select file with IH protocol');
    data_IH = readtable([path_IH filename_IH]);
    data_IH_double = table2array(data_IH);
    tf = data_IH_double(:,1);
    O2_f = data_IH_double(:,2);
    tend = tf(end);
    f0_in = O2_f(1);
end
dt = 0.01;
s = round(tend/dt);
t = zeros(1,s);
dx = 1; dy = 1; % 2-D element dimensions
Vw1 = 1; Vw2 = 1; Vw3 = 1; Vw4 = 1; Vw5 = 1; Vw6 = 1; % (=0 for no liquid present in well and =1 for liquid present)
cells_1 = 1; cells_2 = 1; cells_3 = 1; cells_4 = 1; cells_5 = 1; cells_6 = 1; % (=0 for no cells present in well and =1 for cells present)

% Adjustable system parameters
f0_i = 0.21; % Initial oxygen fraction
Total_volm = 1.5; % Total liquid volume per well (mL)
SC = 0.95; % Salinity correction factor (1 for water and 0.95 for Tyrode's solution)
cO2in = f0_in*P/(R*T_in); % Inlet oxygen concentration

% System initial conditions
c = zeros([14,9,3,s]); % Oxygen concentration in gas phase
c(1:4,1:5,3,1) = f0_i*P/(R*T(2,2,3)); c(1:4,6:9,3,1) = f0_i*P/(R*T(2,6,3)); c(5:9,1:5,3,1) = f0_i*P/(R*T(6,2,3)); c(5:9,6:9,3,1) = f0_i*P/(R*T(6,6,3)); c(10:14,1:5,3,1) = f0_i*P/(R*T(11,2,3)); c(10:14,6:9,3,1) = f0_i*P/(R*T(11,6,3)); c(:,:,2,1) = f0_i*P/(R*T(1,1,2)); c(:,:,1,1) = f0_i*P/(R*T(1,1,1));

% Average oxygen concetration in gas phase above each well
c1_ch = zeros(1,s); c2_ch = zeros(1,s); c3_ch = zeros(1,s); c4_ch = zeros(1,s); c5_ch = zeros(1,s); c6_ch = zeros(1,s);
c1_ch(1) = mean(c(2:4,2:4,3,1),'all'); c2_ch(1) = mean(c(2:4,6:8,3,1),'all'); c3_ch(1) = mean(c(6:8,2:4,3,1),'all'); c4_ch(1) = mean(c(6:8,6:8,3,1),'all'); c5_ch(1) = mean(c(11:13,2:4,3,1),'all'); c6_ch(1) = mean(c(11:13,6:8,3,1),'all');

Pw = zeros([14,9,3,s]); % Water vapor pressure in gas phase
% Average water vapor pressure in gas phase above each well
Pw1_avg = zeros(1,s); Pw2_avg = zeros(1,s); Pw3_avg = zeros(1,s); Pw4_avg = zeros(1,s); Pw5_avg = zeros(1,s); Pw6_avg = zeros(1,s);
Pw1_avg(1) = mean(Pw(2:4,2:4,3,1),'all'); Pw2_avg(1) = mean(Pw(2:4,6:8,3,1),'all'); Pw3_avg(1) = mean(Pw(6:8,2:4,3,1),'all'); Pw4_avg(1) = mean(Pw(6:8,6:8,3,1),'all'); Pw5_avg(1) = mean(Pw(11:13,2:4,3,1),'all'); Pw6_avg(1) = mean(Pw(11:13,6:8,3,1),'all');

% Liquid volume in wells
Vm1 = zeros(1,s); Vm2 = zeros(1,s); Vm3 = zeros(1,s); Vm4 = zeros(1,s); Vm5 = zeros(1,s); Vm6 = zeros(1,s);
Vm1(1) = Total_volm; Vm2(1) = Total_volm; Vm3(1) = Total_volm; Vm4(1) = Total_volm; Vm5(1) = Total_volm; Vm6(1) = Total_volm;

% Initial liquid temperatures for each well
Temp_well1_initial = T_w(2,2,2); Temp_well2_initial = T_w(2,6,2); Temp_well3_initial = T_w(6,2,2); Temp_well4_initial = T_w(6,6,2); Temp_well5_initial = T_w(11,2,2); Temp_well6_initial = T_w(11,6,2);

% Dissolved oxygen concentration in wells (at 20 discrete points for each well)
c1_1 = zeros(1,s); c1_2 = zeros(1,s); c1_3 = zeros(1,s); c1_4 = zeros(1,s); c1_5 = zeros(1,s); c1_6 = zeros(1,s); c1_7 = zeros(1,s); c1_8 = zeros(1,s); c1_9 = zeros(1,s); c1_10 = zeros(1,s); c1_11 = zeros(1,s); c1_12 = zeros(1,s); c1_13 = zeros(1,s); c1_14 = zeros(1,s); c1_15 = zeros(1,s); c1_16 = zeros(1,s); c1_17 = zeros(1,s); c1_18 = zeros(1,s); c1_19 = zeros(1,s); c1_20 = zeros(1,s);
c2_1 = zeros(1,s); c2_2 = zeros(1,s); c2_3 = zeros(1,s); c2_4 = zeros(1,s); c2_5 = zeros(1,s); c2_6 = zeros(1,s); c2_7 = zeros(1,s); c2_8 = zeros(1,s); c2_9 = zeros(1,s); c2_10 = zeros(1,s); c2_11 = zeros(1,s); c2_12 = zeros(1,s); c2_13 = zeros(1,s); c2_14 = zeros(1,s); c2_15 = zeros(1,s); c2_16 = zeros(1,s); c2_17 = zeros(1,s); c2_18 = zeros(1,s); c2_19 = zeros(1,s); c2_20 = zeros(1,s);
c3_1 = zeros(1,s); c3_2 = zeros(1,s); c3_3 = zeros(1,s); c3_4 = zeros(1,s); c3_5 = zeros(1,s); c3_6 = zeros(1,s); c3_7 = zeros(1,s); c3_8 = zeros(1,s); c3_9 = zeros(1,s); c3_10 = zeros(1,s); c3_11 = zeros(1,s); c3_12 = zeros(1,s); c3_13 = zeros(1,s); c3_14 = zeros(1,s); c3_15 = zeros(1,s); c3_16 = zeros(1,s); c3_17 = zeros(1,s); c3_18 = zeros(1,s); c3_19 = zeros(1,s); c3_20 = zeros(1,s);
c4_1 = zeros(1,s); c4_2 = zeros(1,s); c4_3 = zeros(1,s); c4_4 = zeros(1,s); c4_5 = zeros(1,s); c4_6 = zeros(1,s); c4_7 = zeros(1,s); c4_8 = zeros(1,s); c4_9 = zeros(1,s); c4_10 = zeros(1,s); c4_11 = zeros(1,s); c4_12 = zeros(1,s); c4_13 = zeros(1,s); c4_14 = zeros(1,s); c4_15 = zeros(1,s); c4_16 = zeros(1,s); c4_17 = zeros(1,s); c4_18 = zeros(1,s); c4_19 = zeros(1,s); c4_20 = zeros(1,s);
c5_1 = zeros(1,s); c5_2 = zeros(1,s); c5_3 = zeros(1,s); c5_4 = zeros(1,s); c5_5 = zeros(1,s); c5_6 = zeros(1,s); c5_7 = zeros(1,s); c5_8 = zeros(1,s); c5_9 = zeros(1,s); c5_10 = zeros(1,s); c5_11 = zeros(1,s); c5_12 = zeros(1,s); c5_13 = zeros(1,s); c5_14 = zeros(1,s); c5_15 = zeros(1,s); c5_16 = zeros(1,s); c5_17 = zeros(1,s); c5_18 = zeros(1,s); c5_19 = zeros(1,s); c5_20 = zeros(1,s);
c6_1 = zeros(1,s); c6_2 = zeros(1,s); c6_3 = zeros(1,s); c6_4 = zeros(1,s); c6_5 = zeros(1,s); c6_6 = zeros(1,s); c6_7 = zeros(1,s); c6_8 = zeros(1,s); c6_9 = zeros(1,s); c6_10 = zeros(1,s); c6_11 = zeros(1,s); c6_12 = zeros(1,s); c6_13 = zeros(1,s); c6_14 = zeros(1,s); c6_15 = zeros(1,s); c6_16 = zeros(1,s); c6_17 = zeros(1,s); c6_18 = zeros(1,s); c6_19 = zeros(1,s); c6_20 = zeros(1,s);
cells_initial1 = 0.21*(P - interp1(Temp,Psat_range,Temp_well1_initial,'linear','extrap'))*SC*interp1(Temp_sol,Beta_O2_range,Temp_well1_initial,'linear','extrap');
cells_initial2 = 0.21*(P - interp1(Temp,Psat_range,Temp_well2_initial,'linear','extrap'))*SC*interp1(Temp_sol,Beta_O2_range,Temp_well2_initial,'linear','extrap');
cells_initial3 = 0.21*(P - interp1(Temp,Psat_range,Temp_well3_initial,'linear','extrap'))*SC*interp1(Temp_sol,Beta_O2_range,Temp_well3_initial,'linear','extrap');
cells_initial4 = 0.21*(P - interp1(Temp,Psat_range,Temp_well4_initial,'linear','extrap'))*SC*interp1(Temp_sol,Beta_O2_range,Temp_well4_initial,'linear','extrap');
cells_initial5 = 0.21*(P - interp1(Temp,Psat_range,Temp_well5_initial,'linear','extrap'))*SC*interp1(Temp_sol,Beta_O2_range,Temp_well5_initial,'linear','extrap');
cells_initial6 = 0.21*(P - interp1(Temp,Psat_range,Temp_well6_initial,'linear','extrap'))*SC*interp1(Temp_sol,Beta_O2_range,Temp_well6_initial,'linear','extrap');
c1_1(1)= cells_initial1; c1_2(1)= cells_initial1; c1_3(1)= cells_initial1; c1_4(1)= cells_initial1; c1_5(1)= cells_initial1; c1_6(1)= cells_initial1; c1_7(1)= cells_initial1; c1_8(1)= cells_initial1; c1_9(1)= cells_initial1; c1_10(1)= cells_initial1; c1_11(1)= cells_initial1; c1_12(1)= cells_initial1; c1_13(1)= cells_initial1; c1_14(1)= cells_initial1; c1_15(1)= cells_initial1; c1_16(1)= cells_initial1; c1_17(1)= cells_initial1; c1_18(1) = cells_initial1; c1_19(1) = cells_initial1; c1_20(1) = cells_initial1;
c2_1(1)= cells_initial2; c2_2(1)= cells_initial2; c2_3(1)= cells_initial2; c2_4(1)= cells_initial2; c2_5(1)= cells_initial2; c2_6(1)= cells_initial2; c2_7(1)= cells_initial2; c2_8(1)= cells_initial2; c2_9(1)= cells_initial2; c2_10(1)= cells_initial2; c2_11(1)= cells_initial2; c2_12(1)= cells_initial2; c2_13(1)= cells_initial2; c2_14(1)= cells_initial2; c2_15(1)= cells_initial2; c2_16(1)= cells_initial2; c2_17(1)= cells_initial2; c2_18(1) = cells_initial2; c2_19(1) = cells_initial2; c2_20(1) = cells_initial2;
c3_1(1)= cells_initial3; c3_2(1)= cells_initial3; c3_3(1)= cells_initial3; c3_4(1)= cells_initial3; c3_5(1)= cells_initial3; c3_6(1)= cells_initial3; c3_7(1)= cells_initial3; c3_8(1)= cells_initial3; c3_9(1)= cells_initial3; c3_10(1)= cells_initial3; c3_11(1)= cells_initial3; c3_12(1)= cells_initial3; c3_13(1)= cells_initial3; c3_14(1)= cells_initial3; c3_15(1)= cells_initial3; c3_16(1)= cells_initial3; c3_17(1)= cells_initial3; c3_18(1) = cells_initial3; c3_19(1) = cells_initial3; c3_20(1) = cells_initial3;
c4_1(1)= cells_initial4; c4_2(1)= cells_initial4; c4_3(1)= cells_initial4; c4_4(1)= cells_initial4; c4_5(1)= cells_initial4; c4_6(1)= cells_initial4; c4_7(1)= cells_initial4; c4_8(1)= cells_initial4; c4_9(1)= cells_initial4; c4_10(1)= cells_initial4; c4_11(1)= cells_initial4; c4_12(1)= cells_initial4; c4_13(1)= cells_initial4; c4_14(1)= cells_initial4; c4_15(1)= cells_initial4; c4_16(1)= cells_initial4; c4_17(1)= cells_initial4; c4_18(1) = cells_initial4; c4_19(1) = cells_initial4; c4_20(1) = cells_initial4;
c5_1(1)= cells_initial5; c5_2(1)= cells_initial5; c5_3(1)= cells_initial5; c5_4(1)= cells_initial5; c5_5(1)= cells_initial5; c5_6(1)= cells_initial5; c5_7(1)= cells_initial5; c5_8(1)= cells_initial5; c5_9(1)= cells_initial5; c5_10(1)= cells_initial5; c5_11(1)= cells_initial5; c5_12(1)= cells_initial5; c5_13(1)= cells_initial5; c5_14(1)= cells_initial5; c5_15(1)= cells_initial5; c5_16(1)= cells_initial5; c5_17(1)= cells_initial5; c5_18(1) = cells_initial5; c5_19(1) = cells_initial5; c5_20(1) = cells_initial5;
c6_1(1)= cells_initial6; c6_2(1)= cells_initial6; c6_3(1)= cells_initial6; c6_4(1)= cells_initial6; c6_5(1)= cells_initial6; c6_6(1)= cells_initial6; c6_7(1)= cells_initial6; c6_8(1)= cells_initial6; c6_9(1)= cells_initial6; c6_10(1)= cells_initial6; c6_11(1)= cells_initial6; c6_12(1)= cells_initial6; c6_13(1)= cells_initial6; c6_14(1)= cells_initial6; c6_15(1)= cells_initial6; c6_16(1)= cells_initial6; c6_17(1)= cells_initial6; c6_18(1) = cells_initial6; c6_19(1) = cells_initial6; c6_20(1) = cells_initial6;

% Distance between discrete points in well liquid levels
delz1 = zeros(1,s); delz2 = zeros(1,s); delz3 = zeros(1,s); delz4 = zeros(1,s); delz5 = zeros(1,s); delz6 = zeros(1,s);
delz1(1) = (Total_volm/S_w)/19; delz2(1) = (Total_volm/S_w)/19; delz3(1) = (Total_volm/S_w)/19; delz4(1) = (Total_volm/S_w)/19; delz5(1) = (Total_volm/S_w)/19; delz6(1) = (Total_volm/S_w)/19;


%% Solving loop for time-dependent oxygen concentration profile in the chamber
for j=2:round(s)
    t(j) = t(j-1) + dt;

    %% Temperature Data
    % Well 1
    t_m1 = t_m1 + dt; t_curve = t_m1;
    if t_curve>=0 && t_curve<1800
        T(1:4,1:5,3) = 0.0015*t_curve + 20.7;
    elseif t_curve>=1800 && t_curve<4500
        T(1:4,1:5,3) = (17/27000)*t_curve + 334/15;
    elseif t_curve>=4500 && t_curve<7200
        T(1:4,1:5,3) = (1/6750)*t_curve + 733/30;
    else
        T(1:4,1:5,3) = 25.5;
    end

    T_w(2:4,2:4,2) = 1.1*T(2,2,3) - 2.11;

    % Well 2
    if t_curve>=0 && t_curve<1800
        T(1:4,6:9,3) = (7/3000)*t_curve + 104/5;
    elseif t_curve>=1800 && t_curve<4500
        T(1:4,6:9,3) = (7/27000)*t_curve + 368/15;
    elseif t_curve>=4500 && t_curve<7200
        T(1:4,6:9,3) = (1/5400)*t_curve + 373/15;
    else
        T(1:4,6:9,3) = 26.2;
    end

    T_w(2:4,6:8,2) = 1.1*T(2,6,3) - 2.11;

    % Well 3
    if t_curve>=0 && t_curve<1800
        T(5:9,1:5,3) = (83/18000)*t_curve + 103/5;
    elseif t_curve>=1800 && t_curve<4500
        T(5:9,1:5,3) = (17/27000)*t_curve + 833/30;
    elseif t_curve>=4500 && t_curve<7200
        T(5:9,1:5,3) = -(0.2/2700)*t_curve + 83520/2700;
    else
        T(5:9,1:5,3) = 30.4;
    end

    T_w(6:8,2:4,2) = 1.1*T(6,2,3) - 2.11;

    % Well 4
    if t_curve>=0 && t_curve<1800
        T(5:9,6:9,3) = (31/6000)*t_curve + 103/5;
    elseif t_curve>=1800 && t_curve<4500
        T(5:9,6:9,3) = (7/13500)*t_curve + 869/30;
    else
        T(5:9,6:9,3) = 31.3;
    end

    T_w(6:8,6:8,2) = 1.1*T(6,6,3) - 2.11;

    % Well 5
    if t_curve>=0 && t_curve<1800
        T(10:14,1:5,3) = (41/9000)*t_curve + 20.7;
    elseif t_curve>=1800 && t_curve<4500
        T(10:14,1:5,3) = (1/1350)*t_curve + 827/30;
    elseif t_curve>=4500 && t_curve<7200
        T(10:14,1:5,3) = (1/6750)*t_curve + 907/30;
    else
        T(10:14,1:5,3) = 31.3;
    end

    T_w(11:13,2:4,2) = 1.1*T(11,2,3) - 2.11;

    % Well 6
    if t_curve>=0 && t_curve<1800
        T(10:14,6:9,3) = (7/1500)*t_curve + 21;
    elseif t_curve>=1800 && t_curve<4500
        T(10:14,6:9,3) = (17/27000)*t_curve + 424/15;
    elseif t_curve>=4500 && t_curve<7200
        T(10:14,6:9,3) = (1/27000)*t_curve + 464/15;
    else
        T(10:14,6:9,3) = 31.2;
    end

    T_w(11:13,6:8,2) = 1.1*T(11,6,3) - 2.11;

    % Layer 2
    if t_curve>=0 && t_curve<1800
        T(:,:,2) = (7.6/1800)*t_curve + 37800/1800;
    elseif t_curve>=1800 && t_curve<4500
        T(:,:,2) = (2.5/2700)*t_curve + 72720/2700;
    elseif t_curve>=4500 && t_curve<7200
        T(:,:,2) = (0.5/2700)*t_curve + 81720/2700;
    else
        T(:,:,2) = 31.6;
    end

    % Layer 1
    if t_curve>=0 && t_curve<1800
        T(:,:,1) = (10.4/1800)*t_curve + 37440/1800;
    elseif t_curve>=1800 && t_curve<4500
        T(:,:,1) = (1.2/2700)*t_curve + 82080/2700;
    elseif t_curve>=4500 && t_curve<7200
        T(:,:,1) = (1.1/2700)*t_curve + 82530/2700;
    else
        T(:,:,1) = 33.5;
    end

    % Lid temperatures
    % Temperature difference between lid and Layer 3
    T_diff(:,:,1) = -1.3; % Average of experimental data at different positions in the chamber
    if t_curve<1800 % Progression to equilibrium (no temperature difference at experimental t = 0)
        T_diff(:,:,1) = (T_diff(:,:,1)/1800)*t_curve;
    else % Equilibrium reached after 1800s based on measured data
    end
    T_lid = T(:,:,3) + T_diff;
    T = T + 273; T_w(:,:,2) = T_w(:,:,2) + 273; T_lid = T_lid + 273;

    %% Temperature-dependent parameter values
    D1 = D_int(T(2,2,3)); Psat1 = Psat_int(T_w(2,2,2)); Dw1 = Dw_int(T(2,2,3)); DO21 = DO2_int(T_w(2,2,2)); Beta_O21 = Beta_int(T_w(2,2,2)); Psat_cond1 = Psat_int(T_lid(2,2,1));
    D(1:4,1:5,3) = D1;  Psat(2:4,2:4,2) = Psat1;  Dw(1:4,1:5,3) = Dw1;  DO2(2:4,2:4,2) = DO21;  Beta_O2(2:4,2:4,2) = Beta_O21; Psat_cond(1:4,1:5,3) = Psat_cond1;

    D2 = D_int(T(2,6,3)); Psat2 = Psat_int(T_w(2,6,2)); Dw2 = Dw_int(T(2,6,3)); DO22 = DO2_int(T_w(2,6,2)); Beta_O22 = Beta_int(T_w(2,6,2)); Psat_cond2 = Psat_int(T_lid(2,6,1));
    D(1:4,6:9,3) = D2;  Psat(2:4,6:8,2) = Psat2;  Dw(1:4,6:9,3) = Dw2;  DO2(2:4,6:8,2) = DO22;  Beta_O2(2:4,6:8,2) = Beta_O22; Psat_cond(1:4,6:9,3) = Psat_cond2;

    D3 = D_int(T(6,2,3)); Psat3 = Psat_int(T_w(6,2,2)); Dw3 = Dw_int(T(6,2,3)); DO23 = DO2_int(T_w(6,2,2)); Beta_O23 = Beta_int(T_w(6,2,2)); Psat_cond3 = Psat_int(T_lid(6,2,1));
    D(5:9,1:5,3) = D3;  Psat(6:8,2:4,2) = Psat3;  Dw(5:9,1:5,3) = Dw3;  DO2(6:8,2:4,2) = DO23;  Beta_O2(6:8,2:4,2) = Beta_O23; Psat_cond(5:9,1:5,3) = Psat_cond3;

    D4 = D_int(T(6,6,3)); Psat4 = Psat_int(T_w(6,6,2)); Dw4 = Dw_int(T(6,6,3)); DO24 = DO2_int(T_w(6,6,2)); Beta_O24 = Beta_int(T_w(6,6,2)); Psat_cond4 = Psat_int(T_lid(6,6,1));
    D(5:9,6:9,3) = D4;  Psat(6:8,6:8,2) = Psat4;  Dw(5:9,6:9,3) = Dw4;  DO2(6:8,6:8,2) = DO24;  Beta_O2(6:8,6:8,2) = Beta_O24; Psat_cond(5:9,6:9,3) = Psat_cond4;

    D5 = D_int(T(11,2,3)); Psat5 = Psat_int(T_w(11,2,2)); Dw5 = Dw_int(T(11,2,3)); DO25 = DO2_int(T_w(11,2,2)); Beta_O25 = Beta_int(T_w(11,2,2)); Psat_cond5 = Psat_int(T_lid(11,2,1));
    D(10:14,1:5,3) = D5;  Psat(11:13,2:4,2) = Psat5;  Dw(10:14,1:5,3) = Dw5;  DO2(11:13,2:4,2) = DO25;  Beta_O2(11:13,2:4,2) = Beta_O25; Psat_cond(10:14,1:5,3) = Psat_cond5;

    D6 = D_int(T(11,6,3)); Psat6 = Psat_int(T_w(11,6,2)); Dw6 = Dw_int(T(11,6,3)); DO26 = DO2_int(T_w(11,6,2)); Beta_O26 = Beta_int(T_w(11,6,2)); Psat_cond6 = Psat_int(T_lid(11,6,1));
    D(10:14,6:9,3) = D6;  Psat(11:13,6:8,2) = Psat6;  Dw(10:14,6:9,3) = Dw6;  DO2(11:13,6:8,2) = DO26;  Beta_O2(11:13,6:8,2) = Beta_O26; Psat_cond(10:14,6:9,3) = Psat_cond6;

    D1_O2 = DO2(2,2,2); D2_O2 = DO2(2,6,2); D3_O2 = DO2(6,2,2); D4_O2 = DO2(6,6,2); D5_O2 = DO2(11,2,2); D6_O2 = DO2(11,6,2);

    % Layer 1 and 2 diffusion coefficients
    DL1 = D_int(T(1,1,1)); D(:,:,1) = DL1;  DwL1 = Dw_int(T(1,1,1)); Dw(:,:,1) = DwL1;
    DL2 = D_int(T(1,1,2)); D(:,:,2) = DL2;  DwL2 = Dw_int(T(1,1,2)); Dw(:,:,2) = DwL2;

    %% WATER VAPOR PRESSURE AND EVAPORATION MODEL

    % START CODE
    % Layer 2
    z = 2;
    for y=1:14
        if y == 1
            for x=1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dx*dy)*(Pw(y,x,z-1,j-1)-Pw(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + Dw(y,x,z)*(dx*dy)*(Pw(y,x,z+1,j-1)-Pw(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + Pw(y,x,z,j-1);
                elseif x == 2
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + Dw(y,x,z)*(dx*dy)*(Pw(y,x,z-1,j-1)-Pw(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + Dw(y,x,z)*(dx*dy)*(Pw(y,x,z+1,j-1)-Pw(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + Pw(y,x,z,j-1);
                elseif x == 3
                    Pw(y,x,z,j) = 0; % Water vapor pressure in chamber inlet
                elseif x>3 && x<5 || x>5 && x<9
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + Dw(y,x,z)*(dx*dy)*(Pw(y,x,z-1,j-1)-Pw(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + Dw(y,x,z)*(dx*dy)*(Pw(y,x,z+1,j-1)-Pw(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + Pw(y,x,z,j-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dx*dy)*(Pw(y,x,z-1,j-1)-Pw(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + Dw(y,x,z)*(dx*dy)*(Pw(y,x,z+1,j-1)-Pw(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + Pw(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dx*dy)*(Pw(y,x,z-1,j-1)-Pw(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + Dw(y,x,z)*(dx*dy)*(Pw(y,x,z+1,j-1)-Pw(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + Pw(y,x,z,j-1);
                end
            end
        elseif y == 2
            for x=1:9
                if x == 1 || x == 9
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dx*dy)*(Pw(y,x,z-1,j-1)-Pw(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + Dw(y,x,z)*(dx*dy)*(Pw(y,x,z+1,j-1)-Pw(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + Pw(y,x,z,j-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dx*dy)*(Pw(y,x,z+1,j-1)-Pw(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + Pw(y,x,z,j-1);
                else
                    Pw(y,x,z,j) = 0;
                end
            end
        elseif y>2 && y<5 || y>5 && y<9 || y>10 && y<14
            for x=1:9
                if x == 1 || x == 5 || x == 9
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy)*dt/V + Pw(y,x,z,j-1);
                else
                    Pw(y,x,z,j) = 0;
                end
            end
        elseif y == 5
            for x = 1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + Dw(y,x,z)*(dx*(h(y,x,z)))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy)*(dt/V) + Pw(y,x,z,j-1);
                elseif x>1 && x<5 || x>5 && x<9
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x-1,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y-1,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy)*(dt/V) + Pw(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y-1,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy)*(dt/V) + Pw(y,x,z,j-1);
                end
            end
        elseif y == 9
            for x = 1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy)*(dt/V) + Pw(y,x,z,j-1);
                elseif x>1 && x<5 || x>5 && x<9
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x-1,z)*x_y_por_l2(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y-1,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy)*(dt/V) + Pw(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y-1,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l2(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy)*(dt/V) + Pw(y,x,z,j-1);
                end
            end
        elseif y == 10
            for x = 1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 2 || x == 6
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x-1,z)*x_y_por_l2(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y-1,x,z)*(1 - x_y_por_l2(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x-1,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y-1,x,z)*(1 - x_y_por_l2(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 9
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y-1,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy)*(dt/V) + Pw(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x-1,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y-1,x,z)*(1 - x_y_por_l2(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx)*(dt/V) + Pw(y,x,z,j-1);
                end
            end
        else
            for x = 1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx)*(dt/V) + Pw(y,x,z,j-1);
                elseif x>1 && x<5 || x == 6
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x-1,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y-1,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x-1,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y,x+1,z)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx + Q(y,x,z-1)*Pw(y,x,z-1,j-1) - Dw(y,x,z)*(dx*dy)*(Pw(y,x,z,j-1)-Pw(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + Q(y,x,z+1)*Pw(y,x,z+1,j-1) - Dw(y,x,z)*(dx*dy)*(Pw(y,x,z,j-1)-Pw(y,x,z+1,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*Pw(y,x,z,j-1))*(dt/V) + Pw(y,x,z,j-1); % Chamber outlet
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx)*(dt/V) + Pw(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx)*(dt/V) + Pw(y,x,z,j-1);
                end
            end
        end
    end

    % Layer 1
    z = 1;
    for y=1:14
        if y == 1
            for x=1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Dw(y,x,z)*(dx*dy)*(Pw(y,x,z,j-1)-Pw(y,x,z+1,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 2
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Dw(y,x,z)*(dx*dy)*(Pw(y,x,z,j-1)-Pw(y,x,z+1,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 3
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x,z+1,j-1) - Dw(y,x,z)*(dx*dy)*(Pw(y,x,z,j-1)-Pw(y,x,z+1,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx)*(dt/V) + Pw(y,x,z,j-1);
                elseif x>3 && x<9
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Dw(y,x,z)*(dx*dy)*(Pw(y,x,z,j-1)-Pw(y,x,z+1,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx)*(dt/V) + Pw(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Dw(y,x,z)*(dx*dy)*(Pw(y,x,z,j-1)-Pw(y,x,z+1,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy)*(dt/V) + Pw(y,x,z,j-1);
                end
            end
        elseif y == 2
            for x=1:9
                if x == 1 || x == 9
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Dw(y,x,z)*(dx*dy)*(Pw(y,x,z,j-1)-Pw(y,x,z+1,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy)*(dt/V) + Pw(y,x,z,j-1);
                else
                    Pw(y,x,z,j) = 0;
                end
            end
        elseif y>2 && y<14
            for x=1:9
                if x == 1 || x == 9
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy)*(dt/V) + Pw(y,x,z,j-1);
                else
                    Pw(y,x,z,j) = 0;
                end
            end
        else
            for x = 1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx)*(dt/V) + Pw(y,x,z,j-1);
                elseif x>1 && x<7
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x-1,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y,x+1,z)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*dy)*(Pw(y,x,z+1,j-1)-Pw(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx)*(dt/V) + Pw(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx)*(dt/V) + Pw(y,x,z,j-1);
                end
            end
        end
    end

    % Layer 3
    z = 3;
    for y=1:14
        if y == 1
            for x=1:9
                if x < 6
                    k = k1;
                else
                    k = k2;
                end
                cond = k*(dx*dy)*min(0,Psat_cond(y,x,z) - Pw(y,x,z,j-1));
                if x == 1
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Dw(y,x,z)*(dx*dy)*(Pw(y,x,z,j-1)-Pw(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 2
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Dw(y,x,z)*(dx*dy)*(Pw(y,x,z,j-1)-Pw(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 3
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x,z-1,j-1) - Dw(y,x,z)*(dx*dy)*(Pw(y,x,z,j-1)-Pw(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x>3 && x<9
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Dw(y,x,z)*(dx*dy)*(Pw(y,x,z,j-1)-Pw(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Dw(y,x,z)*(dx*dy)*(Pw(y,x,z,j-1)-Pw(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + cond)*(dt/V) + Pw(y,x,z,j-1);
                end
            end
        elseif y == 2
            for x=1:9
                if x < 6
                    k = k1;
                else
                    k = k2;
                end
                cond = k*(dx*dy)*min(0,Psat_cond(y,x,z) - Pw(y,x,z,j-1));
                if x == 1
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y-1,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Dw(y,x,z)*(dx*dy)*(Pw(y,x,z,j-1)-Pw(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 2
                    V = dx*dy*h(y,x,z);
                    if Vw1 == 1
                        if Vm1(j-1)>0
                            kw1 = k1;
                            Evap = kw1*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*(abs(x_n_avg_l3_2)/(abs(x_n_avg_l3_2) + abs(x_p_avg_l3_2)))*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*(abs(x_n_avg_l3_2)/(abs(x_n_avg_l3_2) + abs(x_p_avg_l3_2)))*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                elseif x == 3
                    V = dx*dy*h(y,x,z);
                    if Vw1 == 1
                        if Vm1(j-1)>0
                            kw1 = k1;
                            Evap = kw1*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                elseif x == 4
                    V = dx*dy*h(y,x,z);
                    if Vw1 == 1
                        if Vm1(j-1)>0
                            kw1 = k1;
                            Evap = kw1*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*(abs(x_p_avg_l3_2)/(abs(x_n_avg_l3_2) + abs(x_p_avg_l3_2)))*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*(abs(x_p_avg_l3_2)/(abs(x_n_avg_l3_2) + abs(x_p_avg_l3_2)))*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Dw(y,x,z)*(dx*dy)*(Pw(y,x,z,j-1)-Pw(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x>5 && x<9
                    V = dx*dy*h(y,x,z);
                    if Vw2 == 1
                        if Vm2(j-1)>0
                            kw2 = k2;
                            Evap = kw2*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                else
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y-1,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Dw(y,x,z)*(dx*dy)*(Pw(y,x,z,j-1)-Pw(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + cond)*(dt/V) + Pw(y,x,z,j-1);
                end
            end
        elseif y == 3
            for x=1:9
                if x < 6
                    k = k1;
                else
                    k = k2;
                end
                cond = k*(dx*dy)*min(0,Psat_cond(y,x,z) - Pw(y,x,z,j-1));
                if x == 1
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x>1 && x<5
                    V = dx*dy*h(y,x,z);
                    if Vw1 == 1
                        if Vm1(j-1)>0
                            kw1 = k1;
                            Evap = kw1*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 6
                    V = dx*dy*h(y,x,z);
                    if Vw2 == 1
                        if Vm2(j-1)>0
                            kw2 = k2;
                            Evap = kw2*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    if Vw2 == 1
                        if Vm2(j-1)>0
                            kw2 = k2;
                            Evap = kw2*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    if Vw2 == 1
                        if Vm2(j-1)>0
                            kw2 = k2;
                            Evap = kw2*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                else
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                end
            end
        elseif y == 4
            for x=1:9
                if x < 6
                    k = k1;
                else
                    k = k2;
                end
                cond = k*(dx*dy)*min(0,Psat_cond(y,x,z) - Pw(y,x,z,j-1));
                if x == 1
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x>1 && x<5
                    V = dx*dy*h(y,x,z);
                    if Vw1 == 1
                        if Vm1(j-1)>0
                            kw1 = k1;
                            Evap = kw1*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 6
                    V = dx*dy*h(y,x,z);
                    if Vw2 == 1
                        if Vm2(j-1)>0
                            kw2 = k2;
                            Evap = kw2*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    if Vw2 == 1
                        if Vm2(j-1)>0
                            kw2 = k2;
                            Evap = kw2*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    if Vw2 == 1
                        if Vm2(j-1)>0
                            kw2 = k2;
                            Evap = kw2*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                else
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                end
            end
        elseif y>5 && y<9
            for x=1:9
                if x < 6
                    k = k3;
                else
                    k = k4;
                end
                cond = k*(dx*dy)*min(0,Psat_cond(y,x,z) - Pw(y,x,z,j-1));
                if x == 1
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x>1 && x<5
                    V = dx*dy*h(y,x,z);
                    if Vw3 == 1
                        if Vm3(j-1)>0
                            kw3 = k3;
                            Evap = kw3*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 6
                    V = dx*dy*h(y,x,z);
                    if Vw4 == 1
                        if Vm4(j-1)>0
                            kw4 = k4;
                            Evap = kw4*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    if Vw4 == 1
                        if Vm4(j-1)>0
                            kw4 = k4;
                            Evap = kw4*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    if Vw4 == 1
                        if Vm4(j-1)>0
                            kw4 = k4;
                            Evap = kw4*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                else
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                end
            end
        elseif y>10 && y<14
            for x=1:9
                if x < 6
                    k = k5;
                else
                    k = k6;
                end
                cond = k*(dx*dy)*min(0,Psat_cond(y,x,z) - Pw(y,x,z,j-1));
                if x == 1
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x>1 && x<5
                    V = dx*dy*h(y,x,z);
                    if Vw5 == 1
                        if Vm5(j-1)>0
                            kw5 = k5;
                            Evap = kw5*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 6
                    V = dx*dy*h(y,x,z);
                    if Vw6 == 1
                        if Vm6(j-1)>0
                            kw6 = k6;
                            Evap = kw6*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    if Vw6 == 1
                        if Vm6(j-1)>0
                            kw6 = k6;
                            Evap = kw6*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    if Vw6 == 1
                        if Vm6(j-1)>0
                            kw6 = k6;
                            Evap = kw6*(dx*dy)*(Pw(y,x,z,j-1)-Psat(y,x,2));
                        else
                            Evap = 0;
                        end
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Evap - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    else
                        Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                    end
                else
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                end
            end
        elseif y == 5 || y == 9 || y == 10
            for x=1:9
                if y == 5 || y == 9
                    if x < 6
                        k = k3;
                    else
                        k = k4;
                    end
                else
                    if x < 6
                        k = k5;
                    else
                        k = k6;
                    end
                end
                cond = k*(dx*dy)*min(0,Psat_cond(y,x,z) - Pw(y,x,z,j-1));
                if x == 1
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x>1 && x<7
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y-1,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y+1,x,z,j-1)-Pw(y,x,z,j-1))/dy + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                end
            end
        else
            for x=1:9
                if x < 6
                    k = k5;
                else
                    k = k6;
                end
                cond = k*(dx*dy)*min(0,Psat_cond(y,x,z) - Pw(y,x,z,j-1));
                if x == 1
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x>1 && x<7
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x+1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y-1,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x-1,z)*Pw(y,x-1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x-1,z,j-1))/dx + Q(y,x+1,z)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dx*dy)*(Pw(y,x,z-1,j-1)-Pw(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + cond)*(dt/V) + Pw(y,x,z,j-1);
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy + Q(y,x+1,z)*Pw(y,x+1,z,j-1) - Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y,x+1,z,j-1))/dx - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    Pw(y,x,z,j) = (Q(y,x,z)*Pw(y-1,x,z,j-1) - Dw(y,x,z)*(dx*h(y,x,z))*(Pw(y,x,z,j-1)-Pw(y-1,x,z,j-1))/dy - Q(y,x,z)*Pw(y,x,z,j-1) + Dw(y,x,z)*(dy*h(y,x,z))*(Pw(y,x-1,z,j-1)-Pw(y,x,z,j-1))/dx + cond)*(dt/V) + Pw(y,x,z,j-1);
                end
            end
        end
    end

    % Average water vapor pressure above wells
    Pw1_avg(j) = mean(Pw(2:4,2:4,3,j),'all');
    Pw2_avg(j) = mean(Pw(2:4,6:8,3,j),'all');
    Pw3_avg(j) = mean(Pw(6:8,2:4,3,j),'all');
    Pw4_avg(j) = mean(Pw(6:8,6:8,3,j),'all');
    Pw5_avg(j) = mean(Pw(11:13,2:4,3,j),'all');
    Pw6_avg(j) = mean(Pw(11:13,6:8,3,j),'all');

    % Liquid volume calculation for wells
    if Vw1 == 1
        kc = k1;
        Vm1(j) = -kc*Mw*S_w*dt*(Psat(2,2,2)-Pw1_avg(j-1))/(rho*R*T(2,2,3)) + Vm1(j-1);
    else
        Vm1(j) = 0;
    end

    if Vw2 == 1
        kc = k2;
        Vm2(j) = -kc*Mw*S_w*dt*(Psat(2,6,2)-Pw2_avg(j-1))/(rho*R*T(2,6,3)) + Vm2(j-1);
    else
        Vm2(j) = 0;
    end

    if Vw3 == 1
        kc = k3;
        Vm3(j) = -kc*Mw*S_w*dt*(Psat(6,2,2)-Pw3_avg(j-1))/(rho*R*T(6,2,3)) + Vm3(j-1);
    else
        Vm3(j) = 0;
    end

    if Vw4 == 1
        kc = k4;
        Vm4(j) = -kc*Mw*S_w*dt*(Psat(6,6,2)-Pw4_avg(j-1))/(rho*R*T(6,6,3)) + Vm4(j-1);
    else
        Vm4(j) = 0;
    end

    if Vw5 == 1
        kc = k5;
        Vm5(j) = -kc*Mw*S_w*dt*(Psat(11,2,2)-Pw5_avg(j-1))/(rho*R*T(11,2,3)) + Vm5(j-1);
    else
        Vm5(j) = 0;
    end

    if Vw6 == 1
        kc = k6;
        Vm6(j) = -kc*Mw*S_w*dt*(Psat(11,6,2)-Pw6_avg(j-1))/(rho*R*T(11,6,3)) + Vm6(j-1);
    else
        Vm6(j) = 0;
    end

    % Partition size (distance between 2 discrete points) for well liquid levels
    delz1(j) = (Vm1(j)/S_w)/19; %19 is the number of divisions
    delz2(j) = (Vm2(j)/S_w)/19; %19 is the number of divisions
    delz3(j) = (Vm3(j)/S_w)/19; %19 is the number of divisions
    delz4(j) = (Vm4(j)/S_w)/19; %19 is the number of divisions
    delz5(j) = (Vm5(j)/S_w)/19; %19 is the number of divisions
    delz6(j) = (Vm6(j)/S_w)/19; %19 is the number of divisions
    %END CODE

    %% GAS PHASE OXYGEN MODEL

    %START GAS PHASE OXYGEN CODE
    % Layer 2
    z = 2;
    for y=1:14
        if y == 1
            for x=1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dx*dy)*(c(y,x,z-1,j-1)-c(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + D(y,x,z)*(dx*dy)*(c(y,x,z+1,j-1)-c(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + c(y,x,z,j-1);
                elseif x == 2
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,j-1)-c(y,x,z,j-1))/dx + D(y,x,z)*(dx*dy)*(c(y,x,z-1,j-1)-c(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + D(y,x,z)*(dx*dy)*(c(y,x,z+1,j-1)-c(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + c(y,x,z,j-1);
                elseif x == 3
                    if protocol == 1
                        c(y,x,z,j) = interp1(tf,O2_f*(P/(R*T_in)),t(j),'next','extrap');
                    else
                        c(y,x,z,j) = cO2in;
                    end
                elseif x>3 && x<5 || x>5 && x<9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx + D(y,x,z)*(dx*dy)*(c(y,x,z-1,j-1)-c(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + D(y,x,z)*(dx*dy)*(c(y,x,z+1,j-1)-c(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + c(y,x,z,j-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dx*dy)*(c(y,x,z-1,j-1)-c(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + D(y,x,z)*(dx*dy)*(c(y,x,z+1,j-1)-c(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + c(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dx*dy)*(c(y,x,z-1,j-1)-c(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + D(y,x,z)*(dx*dy)*(c(y,x,z+1,j-1)-c(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + c(y,x,z,j-1);
                end
            end
        elseif y == 2
            for x=1:9
                if x == 1 || x == 9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dx*dy)*(c(y,x,z-1,j-1)-c(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + D(y,x,z)*(dx*dy)*(c(y,x,z+1,j-1)-c(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + c(y,x,z,j-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dx*dy)*(c(y,x,z+1,j-1)-c(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*dt/V + c(y,x,z,j-1);
                else
                    c(y,x,z,j) = 0;
                end
            end
        elseif y>2 && y<5 || y>5 && y<9 || y>10 && y<14
            for x=1:9
                if x == 1 || x == 5 || x == 9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*dt/V + c(y,x,z,j-1);
                else
                    c(y,x,z,j) = 0;
                end
            end
        elseif y == 5
            for x = 1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx + D(y,x,z)*(dx*(h(y,x,z)))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                elseif x>1 && x<5 || x>5 && x<9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x-1,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx + Q(y-1,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy + Q(y,x-1,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                end
            end
        elseif y == 9
            for x = 1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                elseif x>1 && x<5 || x>5 && x<9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx + D(y,x,z)*(dy*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x-1,z)*x_y_por_l2(y,x-1)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx + Q(y-1,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l2(y,x-1)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                end
            end
        elseif y == 10
            for x = 1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                elseif x == 2 || x == 6
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x-1,z)*x_y_por_l2(y,x-1)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx + Q(y-1,x,z)*(1 - x_y_por_l2(y-1,x))*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x-1,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx + Q(y-1,x,z)*(1 - x_y_por_l2(y-1,x))*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                elseif x == 9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy + Q(y,x-1,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x-1,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx + Q(y-1,x,z)*(1 - x_y_por_l2(y-1,x))*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                end
            end
        else
            for x = 1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x>1 && x<5 || x == 6
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x-1,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx + Q(y-1,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x-1,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx + Q(y,x+1,z)*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx + Q(y,x,z-1)*c(y,x,z-1,j-1) - D(y,x,z)*(dx*dy)*(c(y,x,z,j-1)-c(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + Q(y,x,z+1)*c(y,x,z+1,j-1) - D(y,x,z)*(dx*dy)*(c(y,x,z,j-1)-c(y,x,z+1,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*c(y,x,z,j-1))*(dt/V) + c(y,x,z,j-1);
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                end
            end
        end
    end

    % Layer 1
    z = 1;
    for y=1:14
        if y == 1
            for x=1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx - D(y,x,z)*(dx*dy)*(c(y,x,z,j-1)-c(y,x,z+1,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                elseif x == 2
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx - D(y,x,z)*(dx*dy)*(c(y,x,z,j-1)-c(y,x,z+1,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x == 3
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x,z+1,j-1) - D(y,x,z)*(dx*dy)*(c(y,x,z,j-1)-c(y,x,z+1,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,j-1)-c(y,x,z,j-1))/dx + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x>3 && x<9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - D(y,x,z)*(dx*dy)*(c(y,x,z,j-1)-c(y,x,z+1,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - D(y,x,z)*(dx*dy)*(c(y,x,z,j-1)-c(y,x,z+1,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                end
            end
        elseif y == 2
            for x=1:9
                if x == 1 || x == 9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - D(y,x,z)*(dx*dy)*(c(y,x,z,j-1)-c(y,x,z+1,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))) - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                else
                    c(y,x,z,j) = 0;
                end
            end
        elseif y>2 && y<14
            for x=1:9
                if x == 1 || x == 9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                else
                    c(y,x,z,j) = 0;
                end
            end
        else
            for x = 1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x>1 && x<7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x-1,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx + Q(y,x+1,z)*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*dy)*(c(y,x,z+1,j-1)-c(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z+1))))*(dt/V) + c(y,x,z,j-1);
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                end
            end
        end
    end

    % Layer 3
    z = 3;
    for y=1:14
        if y == 1
            for x=1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx - D(y,x,z)*(dx*dy)*(c(y,x,z,j-1)-c(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                elseif x == 2
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx - D(y,x,z)*(dx*dy)*(c(y,x,z,j-1)-c(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x == 3
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x,z-1,j-1) - D(y,x,z)*(dx*dy)*(c(y,x,z,j-1)-c(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,j-1)-c(y,x,z,j-1))/dx + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x>3 && x<9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - D(y,x,z)*(dx*dy)*(c(y,x,z,j-1)-c(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - D(y,x,z)*(dx*dy)*(c(y,x,z,j-1)-c(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                end
            end
        elseif y == 2
            for x=1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - D(y,x,z)*(dx*dy)*(c(y,x,z,j-1)-c(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + Q(y,x+1,z)*x_y_por_l3(y,x+1)*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                elseif x == 2
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*(abs(x_n_avg_l3_2)/(abs(x_n_avg_l3_2) + abs(x_p_avg_l3_2)))*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x == 3
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,j-1)-c(y,x,z,j-1))/dx + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x == 4
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*(abs(x_p_avg_l3_2)/(abs(x_n_avg_l3_2) + abs(x_p_avg_l3_2)))*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x == 5
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - D(y,x,z)*(dx*dy)*(c(y,x,z,j-1)-c(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + Q(y,x-1,z)*x_y_por_l3(y,x-1)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x>5 && x<9
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - D(y,x,z)*(dx*dy)*(c(y,x,z,j-1)-c(y,x,z-1,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))) + Q(y,x-1,z)*x_y_por_l3(y,x-1)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                end
            end
        elseif y == 3
            for x=1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dx*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x>1 && x<7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx + Q(y,x+1,z)*x_y_por_l3(y,x+1)*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                end
            end
        elseif y>3 && y<14
            for x=1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x>1 && x<7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy + Q(y,x-1,z)*x_y_por_l3(y,x-1)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx + Q(y,x+1,z)*x_y_por_l3(y,x+1)*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy)*(dt/V) + c(y,x,z,j-1);
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy + Q(y,x+1,z)*x_y_por_l3(y,x+1)*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y+1,x,z,j-1)-c(y,x,z,j-1))/dy + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                end
            end
        else
            for x=1:9
                if x == 1
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x>1 && x<7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy + Q(y,x-1,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x+1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                elseif x == 7
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy + Q(y,x-1,z)*c(y,x-1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x-1,z,j-1))/dx + Q(y,x+1,z)*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dx*dy)*(c(y,x,z-1,j-1)-c(y,x,z,j-1))/(0.5*(h(y,x,z)+h(y,x,z-1))))*(dt/V) + c(y,x,z,j-1);
                elseif x == 8
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y-1,x,z)*(1-x_y_por_l3(y-1,x))*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy + Q(y,x+1,z)*c(y,x+1,z,j-1) - D(y,x,z)*(dy*h(y,x,z))*(c(y,x,z,j-1)-c(y,x+1,z,j-1))/dx - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                else
                    V = dx*dy*h(y,x,z);
                    c(y,x,z,j) = (Q(y,x,z)*c(y-1,x,z,j-1) - D(y,x,z)*(dx*h(y,x,z))*(c(y,x,z,j-1)-c(y-1,x,z,j-1))/dy - Q(y,x,z)*c(y,x,z,j-1) + D(y,x,z)*(dy*h(y,x,z))*(c(y,x-1,z,j-1)-c(y,x,z,j-1))/dx)*(dt/V) + c(y,x,z,j-1);
                end
            end
        end
    end

    % Average gas-phase oxygen concentration over wells
    c1_ch(j) = mean(c(2:4,2:4,3,j),'all');
    c2_ch(j) = mean(c(2:4,6:8,3,j),'all');
    c3_ch(j) = mean(c(6:8,2:4,3,j),'all');
    c4_ch(j) = mean(c(6:8,6:8,3,j),'all');
    c5_ch(j) = mean(c(11:13,2:4,3,j),'all');
    c6_ch(j) = mean(c(11:13,6:8,3,j),'all');

    % Conversion of average concentrations to partial pressure
    P1_ch(j) = c1_ch(j)*R*T(2,2,3);
    P2_ch(j) = c2_ch(j)*R*T(2,6,3);
    P3_ch(j) = c3_ch(j)*R*T(6,2,3);
    P4_ch(j) = c4_ch(j)*R*T(6,6,3);
    P5_ch(j) = c5_ch(j)*R*T(11,2,3);
    P6_ch(j) = c6_ch(j)*R*T(11,6,3);

    %END GAS PHASE OXYGEN CODE

    %% LIQUID-PHASE DISSOLVED OXYGEN CODE

    % START LIQUID-PHASE OXYGEN CODE
    % Oxygen mass transfer coefficients for each well
    kO2_1 = (D(2,2,3)/Dw(2,2,3))*k1;
    kO2_2 = (D(2,6,3)/Dw(2,6,3))*k2;
    kO2_3 = (D(6,2,3)/Dw(6,2,3))*k3;
    kO2_4 = (D(6,6,3)/Dw(6,6,3))*k4;
    kO2_5 = (D(11,2,3)/Dw(11,2,3))*k5;
    kO2_6 = (D(11,6,3)/Dw(11,6,3))*k6;

    % Well 1
    c1_2(j) = (D1_O2*(c1_3(j-1) - 2*c1_2(j-1) + c1_1(j-1))/(delz1(j-1)^2))*dt + c1_2(j-1);
    c1_1(j) = (kO2_1*c1_ch(j)/D1_O2 + c1_2(j)/delz1(j))/(kO2_1/(SC*Beta_O2(2,2,2)*D1_O2*R*T_w(2,2,2)) + 1/delz1(j));  % Surface boundary condition (finite difference to decompose dc1_1/dz, derivative equivalent to flux)
    c1_3(j) = (D1_O2*(c1_4(j-1) - 2*c1_3(j-1) + c1_2(j-1))/(delz1(j-1)^2))*dt + c1_3(j-1);
    c1_4(j) = (D1_O2*(c1_5(j-1) - 2*c1_4(j-1) + c1_3(j-1))/(delz1(j-1)^2))*dt + c1_4(j-1);
    c1_5(j) = (D1_O2*(c1_6(j-1) - 2*c1_5(j-1) + c1_4(j-1))/(delz1(j-1)^2))*dt + c1_5(j-1);
    c1_6(j) = (D1_O2*(c1_7(j-1) - 2*c1_6(j-1) + c1_5(j-1))/(delz1(j-1)^2))*dt + c1_6(j-1);
    c1_7(j) = (D1_O2*(c1_8(j-1) - 2*c1_7(j-1) + c1_6(j-1))/(delz1(j-1)^2))*dt + c1_7(j-1);
    c1_8(j) = (D1_O2*(c1_9(j-1) - 2*c1_8(j-1) + c1_7(j-1))/(delz1(j-1)^2))*dt + c1_8(j-1);
    c1_9(j) = (D1_O2*(c1_10(j-1) - 2*c1_9(j-1) + c1_8(j-1))/(delz1(j-1)^2))*dt + c1_9(j-1);
    c1_10(j) = (D1_O2*(c1_11(j-1) - 2*c1_10(j-1) + c1_9(j-1))/(delz1(j-1)^2))*dt + c1_10(j-1);
    c1_11(j) = (D1_O2*(c1_12(j-1) - 2*c1_11(j-1) + c1_10(j-1))/(delz1(j-1)^2))*dt + c1_11(j-1);
    c1_12(j) = (D1_O2*(c1_13(j-1) - 2*c1_12(j-1) + c1_11(j-1))/(delz1(j-1)^2))*dt + c1_12(j-1);
    c1_13(j) = (D1_O2*(c1_14(j-1) - 2*c1_13(j-1) + c1_12(j-1))/(delz1(j-1)^2))*dt + c1_13(j-1);
    c1_14(j) = (D1_O2*(c1_15(j-1) - 2*c1_14(j-1) + c1_13(j-1))/(delz1(j-1)^2))*dt + c1_14(j-1);
    c1_15(j) = (D1_O2*(c1_16(j-1) - 2*c1_15(j-1) + c1_14(j-1))/(delz1(j-1)^2))*dt + c1_15(j-1);
    c1_16(j) = (D1_O2*(c1_17(j-1) - 2*c1_16(j-1) + c1_15(j-1))/(delz1(j-1)^2))*dt + c1_16(j-1);
    c1_17(j) = (D1_O2*(c1_18(j-1) - 2*c1_17(j-1) + c1_16(j-1))/(delz1(j-1)^2))*dt + c1_17(j-1);
    c1_18(j) = (D1_O2*(c1_19(j-1) - 2*c1_18(j-1) + c1_17(j-1))/(delz1(j-1)^2))*dt + c1_18(j-1);
    c1_19(j) = (D1_O2*(c1_20(j-1) - 2*c1_19(j-1) + c1_18(j-1))/(delz1(j-1)^2))*dt + c1_19(j-1);
    if cells_1 == 1
        c1_20(j) = (-(ke+O2consump*(delz1(j)/(D1_O2*S_w))-c1_19(j))+sqrt((ke+O2consump*(delz1(j)/(D1_O2*S_w))-c1_19(j))^2+4*c1_19(j)*ke))/2; % Cellular metabolism
    else
        c1_20(j) = c1_19(j); % No flux boundary condition
    end

    % Well 2
    c2_2(j) = (D2_O2*(c2_3(j-1) - 2*c2_2(j-1) + c2_1(j-1))/(delz2(j-1)^2))*dt + c2_2(j-1);
    c2_1(j) = (kO2_2*c2_ch(j)/D2_O2 + c2_2(j)/delz2(j))/(kO2_2/(SC*Beta_O2(2,6,2)*D2_O2*R*T_w(2,6,2)) + 1/delz2(j));
    c2_3(j) = (D2_O2*(c2_4(j-1) - 2*c2_3(j-1) + c2_2(j-1))/(delz2(j-1)^2))*dt + c2_3(j-1);
    c2_4(j) = (D2_O2*(c2_5(j-1) - 2*c2_4(j-1) + c2_3(j-1))/(delz2(j-1)^2))*dt + c2_4(j-1);
    c2_5(j) = (D2_O2*(c2_6(j-1) - 2*c2_5(j-1) + c2_4(j-1))/(delz2(j-1)^2))*dt + c2_5(j-1);
    c2_6(j) = (D2_O2*(c2_7(j-1) - 2*c2_6(j-1) + c2_5(j-1))/(delz2(j-1)^2))*dt + c2_6(j-1);
    c2_7(j) = (D2_O2*(c2_8(j-1) - 2*c2_7(j-1) + c2_6(j-1))/(delz2(j-1)^2))*dt + c2_7(j-1);
    c2_8(j) = (D2_O2*(c2_9(j-1) - 2*c2_8(j-1) + c2_7(j-1))/(delz2(j-1)^2))*dt + c2_8(j-1);
    c2_9(j) = (D2_O2*(c2_10(j-1) - 2*c2_9(j-1) + c2_8(j-1))/(delz2(j-1)^2))*dt + c2_9(j-1);
    c2_10(j) = (D2_O2*(c2_11(j-1) - 2*c2_10(j-1) + c2_9(j-1))/(delz2(j-1)^2))*dt + c2_10(j-1);
    c2_11(j) = (D2_O2*(c2_12(j-1) - 2*c2_11(j-1) + c2_10(j-1))/(delz2(j-1)^2))*dt + c2_11(j-1);
    c2_12(j) = (D2_O2*(c2_13(j-1) - 2*c2_12(j-1) + c2_11(j-1))/(delz2(j-1)^2))*dt + c2_12(j-1);
    c2_13(j) = (D2_O2*(c2_14(j-1) - 2*c2_13(j-1) + c2_12(j-1))/(delz2(j-1)^2))*dt + c2_13(j-1);
    c2_14(j) = (D2_O2*(c2_15(j-1) - 2*c2_14(j-1) + c2_13(j-1))/(delz2(j-1)^2))*dt + c2_14(j-1);
    c2_15(j) = (D2_O2*(c2_16(j-1) - 2*c2_15(j-1) + c2_14(j-1))/(delz2(j-1)^2))*dt + c2_15(j-1);
    c2_16(j) = (D2_O2*(c2_17(j-1) - 2*c2_16(j-1) + c2_15(j-1))/(delz2(j-1)^2))*dt + c2_16(j-1);
    c2_17(j) = (D2_O2*(c2_18(j-1) - 2*c2_17(j-1) + c2_16(j-1))/(delz2(j-1)^2))*dt + c2_17(j-1);
    c2_18(j) = (D2_O2*(c2_19(j-1) - 2*c2_18(j-1) + c2_17(j-1))/(delz2(j-1)^2))*dt + c2_18(j-1);
    c2_19(j) = (D2_O2*(c2_20(j-1) - 2*c2_19(j-1) + c2_18(j-1))/(delz2(j-1)^2))*dt + c2_19(j-1);
    if cells_2 == 1
        c2_20(j) = (-(ke+O2consump*(delz2(j)/(D2_O2*S_w))-c2_19(j))+sqrt((ke+O2consump*(delz2(j)/(D2_O2*S_w))-c2_19(j))^2+4*c2_19(j)*ke))/2; % Cellular metabolism
    else
        c2_20(j) = c2_19(j); % No flux boundary condition
    end

    % Well 3
    c3_2(j) = (D3_O2*(c3_3(j-1) - 2*c3_2(j-1) + c3_1(j-1))/(delz3(j-1)^2))*dt + c3_2(j-1);
    c3_1(j) = (kO2_3*c3_ch(j)/D3_O2 + c3_2(j)/delz3(j))/(kO2_3/(SC*Beta_O2(6,2,2)*D3_O2*R*T_w(6,2,2)) + 1/delz3(j));
    c3_3(j) = (D3_O2*(c3_4(j-1) - 2*c3_3(j-1) + c3_2(j-1))/(delz3(j-1)^2))*dt + c3_3(j-1);
    c3_4(j) = (D3_O2*(c3_5(j-1) - 2*c3_4(j-1) + c3_3(j-1))/(delz3(j-1)^2))*dt + c3_4(j-1);
    c3_5(j) = (D3_O2*(c3_6(j-1) - 2*c3_5(j-1) + c3_4(j-1))/(delz3(j-1)^2))*dt + c3_5(j-1);
    c3_6(j) = (D3_O2*(c3_7(j-1) - 2*c3_6(j-1) + c3_5(j-1))/(delz3(j-1)^2))*dt + c3_6(j-1);
    c3_7(j) = (D3_O2*(c3_8(j-1) - 2*c3_7(j-1) + c3_6(j-1))/(delz3(j-1)^2))*dt + c3_7(j-1);
    c3_8(j) = (D3_O2*(c3_9(j-1) - 2*c3_8(j-1) + c3_7(j-1))/(delz3(j-1)^2))*dt + c3_8(j-1);
    c3_9(j) = (D3_O2*(c3_10(j-1) - 2*c3_9(j-1) + c3_8(j-1))/(delz3(j-1)^2))*dt + c3_9(j-1);
    c3_10(j) = (D3_O2*(c3_11(j-1) - 2*c3_10(j-1) + c3_9(j-1))/(delz3(j-1)^2))*dt + c3_10(j-1);
    c3_11(j) = (D3_O2*(c3_12(j-1) - 2*c3_11(j-1) + c3_10(j-1))/(delz3(j-1)^2))*dt + c3_11(j-1);
    c3_12(j) = (D3_O2*(c3_13(j-1) - 2*c3_12(j-1) + c3_11(j-1))/(delz3(j-1)^2))*dt + c3_12(j-1);
    c3_13(j) = (D3_O2*(c3_14(j-1) - 2*c3_13(j-1) + c3_12(j-1))/(delz3(j-1)^2))*dt + c3_13(j-1);
    c3_14(j) = (D3_O2*(c3_15(j-1) - 2*c3_14(j-1) + c3_13(j-1))/(delz3(j-1)^2))*dt + c3_14(j-1);
    c3_15(j) = (D3_O2*(c3_16(j-1) - 2*c3_15(j-1) + c3_14(j-1))/(delz3(j-1)^2))*dt + c3_15(j-1);
    c3_16(j) = (D3_O2*(c3_17(j-1) - 2*c3_16(j-1) + c3_15(j-1))/(delz3(j-1)^2))*dt + c3_16(j-1);
    c3_17(j) = (D3_O2*(c3_18(j-1) - 2*c3_17(j-1) + c3_16(j-1))/(delz3(j-1)^2))*dt + c3_17(j-1);
    c3_18(j) = (D3_O2*(c3_19(j-1) - 2*c3_18(j-1) + c3_17(j-1))/(delz3(j-1)^2))*dt + c3_18(j-1);
    c3_19(j) = (D3_O2*(c3_20(j-1) - 2*c3_19(j-1) + c3_18(j-1))/(delz3(j-1)^2))*dt + c3_19(j-1);
    if cells_3 == 1
        c3_20(j) = (-(ke+O2consump*(delz3(j)/(D3_O2*S_w))-c3_19(j))+sqrt((ke+O2consump*(delz3(j)/(D3_O2*S_w))-c3_19(j))^2+4*c3_19(j)*ke))/2; % Cellular metabolism
    else
        c3_20(j) = c3_19(j); % No flux boundary condition
    end

    % Well 4
    c4_2(j) = (D4_O2*(c4_3(j-1) - 2*c4_2(j-1) + c4_1(j-1))/(delz4(j-1)^2))*dt + c4_2(j-1);
    c4_1(j) = (kO2_4*c4_ch(j)/D4_O2 + c4_2(j)/delz4(j))/(kO2_4/(SC*Beta_O2(6,6,2)*D4_O2*R*T_w(6,6,2)) + 1/delz4(j));
    c4_3(j) = (D4_O2*(c4_4(j-1) - 2*c4_3(j-1) + c4_2(j-1))/(delz4(j-1)^2))*dt + c4_3(j-1);
    c4_4(j) = (D4_O2*(c4_5(j-1) - 2*c4_4(j-1) + c4_3(j-1))/(delz4(j-1)^2))*dt + c4_4(j-1);
    c4_5(j) = (D4_O2*(c4_6(j-1) - 2*c4_5(j-1) + c4_4(j-1))/(delz4(j-1)^2))*dt + c4_5(j-1);
    c4_6(j) = (D4_O2*(c4_7(j-1) - 2*c4_6(j-1) + c4_5(j-1))/(delz4(j-1)^2))*dt + c4_6(j-1);
    c4_7(j) = (D4_O2*(c4_8(j-1) - 2*c4_7(j-1) + c4_6(j-1))/(delz4(j-1)^2))*dt + c4_7(j-1);
    c4_8(j) = (D4_O2*(c4_9(j-1) - 2*c4_8(j-1) + c4_7(j-1))/(delz4(j-1)^2))*dt + c4_8(j-1);
    c4_9(j) = (D4_O2*(c4_10(j-1) - 2*c4_9(j-1) + c4_8(j-1))/(delz4(j-1)^2))*dt + c4_9(j-1);
    c4_10(j) = (D4_O2*(c4_11(j-1) - 2*c4_10(j-1) + c4_9(j-1))/(delz4(j-1)^2))*dt + c4_10(j-1);
    c4_11(j) = (D4_O2*(c4_12(j-1) - 2*c4_11(j-1) + c4_10(j-1))/(delz4(j-1)^2))*dt + c4_11(j-1);
    c4_12(j) = (D4_O2*(c4_13(j-1) - 2*c4_12(j-1) + c4_11(j-1))/(delz4(j-1)^2))*dt + c4_12(j-1);
    c4_13(j) = (D4_O2*(c4_14(j-1) - 2*c4_13(j-1) + c4_12(j-1))/(delz4(j-1)^2))*dt + c4_13(j-1);
    c4_14(j) = (D4_O2*(c4_15(j-1) - 2*c4_14(j-1) + c4_13(j-1))/(delz4(j-1)^2))*dt + c4_14(j-1);
    c4_15(j) = (D4_O2*(c4_16(j-1) - 2*c4_15(j-1) + c4_14(j-1))/(delz4(j-1)^2))*dt + c4_15(j-1);
    c4_16(j) = (D4_O2*(c4_17(j-1) - 2*c4_16(j-1) + c4_15(j-1))/(delz4(j-1)^2))*dt + c4_16(j-1);
    c4_17(j) = (D4_O2*(c4_18(j-1) - 2*c4_17(j-1) + c4_16(j-1))/(delz4(j-1)^2))*dt + c4_17(j-1);
    c4_18(j) = (D4_O2*(c4_19(j-1) - 2*c4_18(j-1) + c4_17(j-1))/(delz4(j-1)^2))*dt + c4_18(j-1);
    c4_19(j) = (D4_O2*(c4_20(j-1) - 2*c4_19(j-1) + c4_18(j-1))/(delz4(j-1)^2))*dt + c4_19(j-1);
    if cells_4 == 1
        c4_20(j) = (-(ke+O2consump*(delz4(j)/(D4_O2*S_w))-c4_19(j))+sqrt((ke+O2consump*(delz4(j)/(D4_O2*S_w))-c4_19(j))^2+4*c4_19(j)*ke))/2; % Cellular metabolism
    else
        c4_20(j) = c4_19(j); % No flux boundary condition
    end

    % Well 5
    c5_2(j) = (D5_O2*(c5_3(j-1) - 2*c5_2(j-1) + c5_1(j-1))/(delz5(j-1)^2))*dt + c5_2(j-1);
    c5_1(j) = (kO2_5*c5_ch(j)/D5_O2 + c5_2(j)/delz5(j))/(kO2_5/(SC*Beta_O2(11,2,2)*D5_O2*R*T_w(11,2,2)) + 1/delz5(j));
    c5_3(j) = (D5_O2*(c5_4(j-1) - 2*c5_3(j-1) + c5_2(j-1))/(delz5(j-1)^2))*dt + c5_3(j-1);
    c5_4(j) = (D5_O2*(c5_5(j-1) - 2*c5_4(j-1) + c5_3(j-1))/(delz5(j-1)^2))*dt + c5_4(j-1);
    c5_5(j) = (D5_O2*(c5_6(j-1) - 2*c5_5(j-1) + c5_4(j-1))/(delz5(j-1)^2))*dt + c5_5(j-1);
    c5_6(j) = (D5_O2*(c5_7(j-1) - 2*c5_6(j-1) + c5_5(j-1))/(delz5(j-1)^2))*dt + c5_6(j-1);
    c5_7(j) = (D5_O2*(c5_8(j-1) - 2*c5_7(j-1) + c5_6(j-1))/(delz5(j-1)^2))*dt + c5_7(j-1);
    c5_8(j) = (D5_O2*(c5_9(j-1) - 2*c5_8(j-1) + c5_7(j-1))/(delz5(j-1)^2))*dt + c5_8(j-1);
    c5_9(j) = (D5_O2*(c5_10(j-1) - 2*c5_9(j-1) + c5_8(j-1))/(delz5(j-1)^2))*dt + c5_9(j-1);
    c5_10(j) = (D5_O2*(c5_11(j-1) - 2*c5_10(j-1) + c5_9(j-1))/(delz5(j-1)^2))*dt + c5_10(j-1);
    c5_11(j) = (D5_O2*(c5_12(j-1) - 2*c5_11(j-1) + c5_10(j-1))/(delz5(j-1)^2))*dt + c5_11(j-1);
    c5_12(j) = (D5_O2*(c5_13(j-1) - 2*c5_12(j-1) + c5_11(j-1))/(delz5(j-1)^2))*dt + c5_12(j-1);
    c5_13(j) = (D5_O2*(c5_14(j-1) - 2*c5_13(j-1) + c5_12(j-1))/(delz5(j-1)^2))*dt + c5_13(j-1);
    c5_14(j) = (D5_O2*(c5_15(j-1) - 2*c5_14(j-1) + c5_13(j-1))/(delz5(j-1)^2))*dt + c5_14(j-1);
    c5_15(j) = (D5_O2*(c5_16(j-1) - 2*c5_15(j-1) + c5_14(j-1))/(delz5(j-1)^2))*dt + c5_15(j-1);
    c5_16(j) = (D5_O2*(c5_17(j-1) - 2*c5_16(j-1) + c5_15(j-1))/(delz5(j-1)^2))*dt + c5_16(j-1);
    c5_17(j) = (D5_O2*(c5_18(j-1) - 2*c5_17(j-1) + c5_16(j-1))/(delz5(j-1)^2))*dt + c5_17(j-1);
    c5_18(j) = (D5_O2*(c5_19(j-1) - 2*c5_18(j-1) + c5_17(j-1))/(delz5(j-1)^2))*dt + c5_18(j-1);
    c5_19(j) = (D5_O2*(c5_20(j-1) - 2*c5_19(j-1) + c5_18(j-1))/(delz5(j-1)^2))*dt + c5_19(j-1);
    if cells_5 == 1
        c5_20(j) = (-(ke+O2consump*(delz5(j)/(D5_O2*S_w))-c5_19(j))+sqrt((ke+O2consump*(delz5(j)/(D5_O2*S_w))-c5_19(j))^2+4*c5_19(j)*ke))/2; % Cellular metabolism
    else
        c5_20(j) = c5_19(j); % No flux boundary condition
    end

    % Well 6
    c6_2(j) = (D6_O2*(c6_3(j-1) - 2*c6_2(j-1) + c6_1(j-1))/(delz6(j-1)^2))*dt + c6_2(j-1);
    c6_1(j) = (kO2_6*c6_ch(j)/D6_O2 + c6_2(j)/delz6(j))/(kO2_6/(SC*Beta_O2(11,6,2)*D6_O2*R*T_w(11,6,2)) + 1/delz6(j));
    c6_3(j) = (D6_O2*(c6_4(j-1) - 2*c6_3(j-1) + c6_2(j-1))/(delz6(j-1)^2))*dt + c6_3(j-1);
    c6_4(j) = (D6_O2*(c6_5(j-1) - 2*c6_4(j-1) + c6_3(j-1))/(delz6(j-1)^2))*dt + c6_4(j-1);
    c6_5(j) = (D6_O2*(c6_6(j-1) - 2*c6_5(j-1) + c6_4(j-1))/(delz6(j-1)^2))*dt + c6_5(j-1);
    c6_6(j) = (D6_O2*(c6_7(j-1) - 2*c6_6(j-1) + c6_5(j-1))/(delz6(j-1)^2))*dt + c6_6(j-1);
    c6_7(j) = (D6_O2*(c6_8(j-1) - 2*c6_7(j-1) + c6_6(j-1))/(delz6(j-1)^2))*dt + c6_7(j-1);
    c6_8(j) = (D6_O2*(c6_9(j-1) - 2*c6_8(j-1) + c6_7(j-1))/(delz6(j-1)^2))*dt + c6_8(j-1);
    c6_9(j) = (D6_O2*(c6_10(j-1) - 2*c6_9(j-1) + c6_8(j-1))/(delz6(j-1)^2))*dt + c6_9(j-1);
    c6_10(j) = (D6_O2*(c6_11(j-1) - 2*c6_10(j-1) + c6_9(j-1))/(delz6(j-1)^2))*dt + c6_10(j-1);
    c6_11(j) = (D6_O2*(c6_12(j-1) - 2*c6_11(j-1) + c6_10(j-1))/(delz6(j-1)^2))*dt + c6_11(j-1);
    c6_12(j) = (D6_O2*(c6_13(j-1) - 2*c6_12(j-1) + c6_11(j-1))/(delz6(j-1)^2))*dt + c6_12(j-1);
    c6_13(j) = (D6_O2*(c6_14(j-1) - 2*c6_13(j-1) + c6_12(j-1))/(delz6(j-1)^2))*dt + c6_13(j-1);
    c6_14(j) = (D6_O2*(c6_15(j-1) - 2*c6_14(j-1) + c6_13(j-1))/(delz6(j-1)^2))*dt + c6_14(j-1);
    c6_15(j) = (D6_O2*(c6_16(j-1) - 2*c6_15(j-1) + c6_14(j-1))/(delz6(j-1)^2))*dt + c6_15(j-1);
    c6_16(j) = (D6_O2*(c6_17(j-1) - 2*c6_16(j-1) + c6_15(j-1))/(delz6(j-1)^2))*dt + c6_16(j-1);
    c6_17(j) = (D6_O2*(c6_18(j-1) - 2*c6_17(j-1) + c6_16(j-1))/(delz6(j-1)^2))*dt + c6_17(j-1);
    c6_18(j) = (D6_O2*(c6_19(j-1) - 2*c6_18(j-1) + c6_17(j-1))/(delz6(j-1)^2))*dt + c6_18(j-1);
    c6_19(j) = (D6_O2*(c6_20(j-1) - 2*c6_19(j-1) + c6_18(j-1))/(delz6(j-1)^2))*dt + c6_19(j-1);
    if cells_6 == 1
        c6_20(j) = (-(ke+O2consump*(delz6(j)/(D6_O2*S_w))-c6_19(j))+sqrt((ke+O2consump*(delz6(j)/(D6_O2*S_w))-c6_19(j))^2+4*c6_19(j)*ke))/2; % Cellular metabolism
    else
        c6_20(j) = c6_19(j); % No flux boundary condition
    end

end
