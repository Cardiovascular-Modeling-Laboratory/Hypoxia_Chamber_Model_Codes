%% Load BHM results and extract and filter data
% Load model results of posterior simulations from "Gas_and_liquid_phase_calibration.ipynb'
[y_post_name, y_post_path] = uigetfile('*.csv', 'Select file with model results of posterior simulations (titled y_out_post_liquid_final.csv)'); 
post_O2 = readtable([y_post_path y_post_name]);
post_O2_data = table2array(post_O2);

% Extract data 
V1_s = post_O2_data(:,1); % Volume of water remainining in Well 1
V2_s = post_O2_data(:,2); % Volume of water remainining in Well 2
V3_s = post_O2_data(:,3); % Volume of water remainining in Well 3
V4_s = post_O2_data(:,4); % Volume of water remainining in Well 4
V5_s = post_O2_data(:,5); % Volume of water remainining in Well 5
V6_s = post_O2_data(:,6); % Volume of water remainining in Well 6
c3_5_s = post_O2_data(:,7); % Dissolved oxygen concentration in Well 3 at assumed probe position after 5 min. of step decrease in oxygen to 0
c3_10_s = post_O2_data(:,8); % Dissolved oxygen concentration in Well 3 at assumed probe position after 10 min. of step decrease in oxygen to 0
c3_15_s = post_O2_data(:,9); % Dissolved oxygen concentration in Well 3 at assumed probe position after 15 min. of step decrease in oxygen to 0

% Filter data (Determine parameter sets which result in the minimum, maximum, and mean for each output variable)
% Initialize
max_p = zeros(1,9);
min_p = zeros(1,9);
mean_p = zeros(1,9);

% Volume of water remaining in Well 1
max_p(1) = find(V1_s == max(V1_s));
min_p(1) = find(V1_s == min(V1_s));
error = abs(V1_s-mean(V1_s));
[~, ind] = mink(error,1);
mean_p(1) = ind;

% Volume of water remaining in Well 2
max_p(2) = find(V2_s == max(V2_s));
min_p(2) = find(V2_s == min(V2_s));
error = abs(V2_s-mean(V2_s));
[~, ind] = mink(error,1);
mean_p(2) = ind;

% Volume of water remaining in Well 3
max_p(3) = find(V3_s == max(V3_s));
min_p(3) = find(V3_s == min(V3_s));
error = abs(V3_s-mean(V3_s));
[~, ind] = mink(error,1);
mean_p(3) = ind;

% Volume of water remaining in Well 4
max_p(4) = find(V4_s == max(V4_s));
min_p(4) = find(V4_s == min(V4_s));
error = abs(V4_s-mean(V4_s));
[~, ind] = mink(error,1);
mean_p(4) = ind;

% Volume of water remaining in Well 5
max_p(5) = find(V5_s == max(V5_s));
min_p(5) = find(V5_s == min(V5_s));
error = abs(V5_s-mean(V5_s));
[~, ind] = mink(error,1);
mean_p(5) = ind;

% Volume of water remaining in Well 6
max_p(6) = find(V6_s == max(V6_s));
min_p(6) = find(V6_s == min(V6_s));
error = abs(V6_s-mean(V6_s));
[~, ind] = mink(error,1);
mean_p(6) = ind;

% Oxygen concentration in Well 3 after 5 min. of step decrease
max_p(7) = find(c3_5_s == max(c3_5_s));
min_p(7) = find(c3_5_s == min(c3_5_s));
error = abs(c3_5_s-mean(c3_5_s));
[~, ind] = mink(error,1);
mean_p(7) = ind;

% Oxygen concentration in Well 3 after 10 min. of step decrease
max_p(8) = find(c3_10_s == max(c3_10_s));
min_p(8) = find(c3_10_s == min(c3_10_s));
error = abs(c3_10_s-mean(c3_10_s));
[~, ind] = mink(error,1);
mean_p(8) = ind;

% Oxygen concentration in Well 3 after 15 min. of step decrease
max_p(9) = find(c3_15_s == max(c3_15_s));
min_p(9) = find(c3_15_s == min(c3_15_s));
error = abs(c3_15_s-mean(c3_15_s));
[~, ind] = mink(error,1);
mean_p(9) = ind;

pars_run_total = [max_p, min_p, mean_p]; 
pars_run = unique(pars_run_total); % Eliminate any repeating parameter sets
[x_post_name, x_post_path] = uigetfile('*.csv', 'Select file with parameter sets used for posterior simulations (titled x_pars_post_liquid_final.csv)'); 
x_pars = readtable([x_post_path x_post_name]);
x_pars_v = table2array(x_pars);
x_pars_filtered = x_pars_v(pars_run,:); % Columns are (T1, k1, k2, k3, k4, k5, k6, Inlet_l2, Inlet_l3, position, td)

