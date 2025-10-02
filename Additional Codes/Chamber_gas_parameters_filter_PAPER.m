%% Load BHM results and extract and filter data
% Load model results of posterior simulations from "Gas_phase_calibration.ipynb'
[y_post_name, y_post_path] = uigetfile('*.csv', 'Select file with model results of posterior simulations (titled y_out_post_gas_final.csv)'); 
post_O2 = readtable([y_post_path y_post_name]);
post_O2_data = table2array(post_O2);
% Extract data 
TR_10_s = post_O2_data(:,1); % Oxygen partial pressure at top right position after 10 s. of step increase in oxygen fraction to 0.21
TR_60_s = post_O2_data(:,2); % Oxygen partial pressure at top right position after 60 s. of step increase in oxygen fraction to 0.21
TR_90_s = post_O2_data(:,3); % Oxygen partial pressure at top right position after 90 s. of step increase in oxygen fraction to 0.21
BL_30_s = post_O2_data(:,4); % Oxygen partial pressure at bottom left position after 30 s. of step increase in oxygen fraction to 0.21
BL_60_s = post_O2_data(:,5); % Oxygen partial pressure at bottom left position after 60 s. of step increase in oxygen fraction to 0.21
BL_90_s = post_O2_data(:,6); % Oxygen partial pressure at bottom left position after 90 s. of step increase in oxygen fraction to 0.21

% Filter data (Determine parameter sets which result in the minimum, maximum, and mean for each output variable)
% Initialize
max_p = zeros(1,6);
min_p = zeros(1,6);
mean_p = zeros(1,6);

% Oxygen partial pressure at top right postion 10 s. after step increase
max_p(1) = find(TR_10_s == max(TR_10_s));
min_p(1) = find(TR_10_s == min(TR_10_s));
error = abs(TR_10_s-mean(TR_10_s));
[~, ind] = mink(error,1);
mean_p(1) = ind;

% Oxygen partial pressure at top right postion 60 s. after step increase
max_TR = find(TR_60_s == max(TR_60_s));
max_p(2) = max_TR(1);
min_p(2) = find(TR_60_s == min(TR_60_s));
error = abs(TR_60_s-mean(TR_60_s));
[~, ind] = mink(error,1);
mean_p(2) = ind;

% Oxygen partial pressure at top right postion 90 s. after step increase
max_p(3) = find(TR_90_s == max(TR_90_s));
min_p(3) = find(TR_90_s == min(TR_90_s));
error = abs(TR_90_s-mean(TR_90_s));
[~, ind] = mink(error,1);
mean_p(3) = ind;

% Oxygen partial pressure at bottom left postion 30 s. after step increase
max_p(4) = find(BL_30_s == max(BL_30_s));
min_p(4) = find(BL_30_s == min(BL_30_s));
error = abs(BL_30_s-mean(BL_30_s));
[~, ind] = mink(error,1);
mean_p(4) = ind;

% Oxygen partial pressure at bottom left postion 60 s. after step increase
max_p(5) = find(BL_60_s == max(BL_60_s));
min_p(5) = find(BL_60_s == min(BL_60_s));
error = abs(BL_60_s-mean(BL_60_s));
[~, ind] = mink(error,1);
mean_p(5) = ind;

% Oxygen partial pressure at bottom left postion 90 s. after step increase
max_p(6) = find(BL_90_s == max(BL_90_s));
min_p(6) = find(BL_90_s == min(BL_90_s));
error = abs(BL_90_s-mean(BL_90_s));
[~, ind] = mink(error,1);
mean_p(6) = ind;

pars_run_total = [max_p, min_p, mean_p];
pars_run = unique(pars_run_total); % Eliminate any repeating parameter sets
[x_post_name, x_post_path] = uigetfile('*.csv', 'Select file with parameter sets used for posterior simulations (titled x_pars_post_gas_final.csv)'); 
x_pars = readtable([x_post_path x_post_name]);
x_pars_v = table2array(x_pars);
x_pars_filtered = x_pars_v(pars_run,:); % Columns are (T1, td, Inlet_l2, Inlet_l3)

