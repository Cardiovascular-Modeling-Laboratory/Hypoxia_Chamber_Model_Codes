%% Load and segment target oxygen pattern
% Load .mat file with target oxygen pattern
[filename_cycle, cycle_path] = ...
    uigetfile('*.mat', 'Select the .mat file with the target oxygen fluctuation pattern',[]);
file_path = [cycle_path filename_cycle];
load(file_path)

% Segment input pattern
diffc_t = diff(c_t)/(t1(2)-t1(1)); % Pattern derivative
indices = find([0 diff(sign(diffc_t))]~=0); % Identify inflection points
indices = [1 indices length(c_t)];

% Identify local extrema for each segment
target_c(1) = c_t(indices(1));
for i=2:length(indices)-1
    sign_i = diffc_t(indices(i)+10); % Determine sign of derivative after inflection point
    % Estimate oxygen concentration at inflection point from input pattern
    if sign_i>0
        target_c(i) = min(CsaO2d_target(indices(i-1):indices(i+1)));
    else
        target_c(i) = max(CsaO2d_target(indices(i-1):indices(i+1)));
    end
end
target_c(2:end) = target_c(2:end)*1e9;
target_c = target_c(2:end);

%% Fit models
% Select file with GPE training data
[filename_GPE, GPE_path] = ...
    uigetfile('*.xlsx', 'Select the Excel file with the GPE training data',[]);
file_path_GPE = [GPE_path filename_GPE];
sheet = 'Test_0_75mL';

% Remnant increase model
range_incx = 'A3:D34';
range_incy = 'F3:F34';
x_inc = readtable(file_path_GPE, 'Sheet', sheet, 'Range', range_incx);
x_inc = table2array(x_inc);
y_inc = readtable(file_path_GPE, 'Sheet', sheet, 'Range', range_incy);
y_inc = table2array(y_inc);
model_inc = fitrgp(x_inc, y_inc, 'KernelFunction', 'squaredexponential', 'Sigma', 0.1, 'Standardize', true); 

% Remnant decrease model
range_decx = 'I3:L34';
range_decy = 'N3:N34';
x_dec = readtable(file_path_GPE, 'Sheet', sheet, 'Range', range_decx);
x_dec = table2array(x_dec);
y_dec = readtable(file_path_GPE, 'Sheet', sheet, 'Range', range_decy);
y_dec = table2array(y_dec);
model_dec = fitrgp(x_dec, y_dec, 'KernelFunction', 'squaredexponential', 'Sigma', 0.1, 'Standardize', true);

%% Prepare parameter set for optimization run (FIRST RUN CHAMBER_LIQUID_PARAMETERS_FILTER_PAPER.M)
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

%% Temperature Inlet
% Well 1
if T1>=20.7 && T1<23.4
    t_m1 = (2000/3)*T1 - 13800;
elseif T1>=23.4 && T1<25.1
    t_m1 = (27000/17)*T1 - 601200/17;
elseif T1>=25.1 && T1<25.5
    t_m1 = 6750*T1 - 164925;
else
    t_m1 = 7200;
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

T_in = L2 + 273; % Estimate chamber inlet for the run (assumed to be constant)

%% Protocol Bounds and Optimization Thresholds
% Inlet oxygen fraction bounds for decreasing oxygen segment
ll_dec = 0;
ul_dec = 0.1;
mid_dec = 0.5*(ll_dec + ul_dec);

% Inlet oxygen fraction bounds for increasing oxygen segment
ll_inc = 0.25;
ul_inc = 0.5;
mid_inc = 0.5*(ll_inc + ul_inc);

% Run duration bounds (in seconds)
ll_time = 4;
ul_time = 60;

% Optimization thresholds
O2_threshold = 2; % uM
time_threshold = 10; % seconds
range_threshold = 5; % For GPE predictions (uM)
O2f_step = 0.01;

%% Optimization Code
% Loop through each segment of the target pattern
for i = 1:length(target_c)-1
    % First iteration based on previous model testing
    if i == 1
        sol_f(i) = 0.1;
        sol_dur(i) = 600;
    elseif i == 2
        % Optimization targets for segment
        target_t = t1(indices(i+1)) - t1(indices(i)); % Target time for segment
        target_prev = target_c(i-1); % Target oxygen concentration at beginning of segment
        target = target_c(i); % Target oxygen concentration at end of segment

        % Run model with results of previous iteration to obtain the initial condition for current iteration
        f0 = sol_f(i-1); % Inlet oxygen fraction from previous iteration
        tdur = sol_dur(i-1); % Run duration from previous iteration
        sol_prev = Chamber_initial_condition(T1,k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,f0,tdur);
        ic = sol_prev; % Initial condition of chamber stored in structure

        % Run current iteration
        clear objfun options
        % Determine whether segment will have increasing or decreasing oxygen concentration
        if target_c(i) > target_c(i-1) % Increasing oxygen concentration
            O2_dir = 0;
        elseif target_c(i) < target_c(i-1) % Decreasing oxygen concentration
            O2_dir = 1;
        else
        end

        % Adjust thresholds based on direction of oxygen profile for segment (time threshold is the same for both)
        if O2_dir == 0 % Increasing oxygen concentration
            lb = [ll_inc, ll_time];
            ub = [ul_inc, ul_time];
        else % Decreasing oxygen concentration
            lb = [ll_dec, ll_time];
            ub = [ul_dec, ul_time];
        end

        % Create optimization function and run Pareto front solver (gamultiobj)
        toggle = 1; % Define function output to be objective functions (toggle = 2 is to have function output be the chamber state for initial condition)
        objfun = @(fdur)opt_function_secondary(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,fdur(1),fdur(2),target_prev,target,target_t,toggle,ic,O2_dir,T_in);
        options = optimoptions('gamultiobj', 'Display','iter','MaxGenerations', 10,'PopulationSize',100,'ParetoFraction',0.5,'UseParallel', true);
        [sol, fval] = gamultiobj(objfun, 2, [], [], [], [], lb, ub, [], options); % Run optimization
        test_f = sol(:,1); % Solutions for inlet oxygen fraction
        test_dur = sol(:,2); % Solutions for run duration

        % Results of objective functions
        time_val = fval(:,1);
        target_val = fval(:,2);
        target_prev_val = fval(:,3);

        % Find solution set from optimization results
        psol_ind = find(target_prev_val < O2_threshold); % Filter possible solution sets based on oxygen concentration threshold
        time_psol = time_val(psol_ind);
        time_psol_sort = unique(sort(time_psol)); % Eliminate repeating solution sets
        sol_row = []; % Initialize to enter while loop
        counter_loop = 1;

        while double(isempty(sol_row)) == 1 % Run through possible solution sets until defined criteria is met
            test_row = find(time_val == time_psol_sort(counter_loop));
            test_row = test_row(1);
            test_x_1 = [ic.c3_20_ic*1e9, test_f(test_row), test_dur(test_row), ll_dec];
            pred_1 = predict(model_inc,test_x_1); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
            test_x_2 = [ic.c3_20_ic*1e9, test_f(test_row), test_dur(test_row), ul_dec];
            pred_2 = predict(model_inc,test_x_2); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment

            if target_val(test_row) > pred_1 - O2_threshold && target_val(test_row) < pred_2 + O2_threshold % Target remnant increase falls within GPE predicted range with uncertainty
                test_x_3 = [ic.c3_20_ic*1e9, test_f(test_row), test_dur(test_row), mid_dec];
                pred_3 = predict(model_inc,test_x_3); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                pred_range = [pred_1, pred_3, pred_2]; % GPE prediction range
                [pred_sorted, idx] = sort(pred_range);
                test_vals = [ll_dec, mid_dec, ul_dec];
                f_sorted = test_vals(idx);
                if length(unique(round(pred_range,3))) < length(pred_range) || range(pred_range) < range_threshold % There does not seem to be significant influence of inlet oxygen fraction for next segment
                    sol_f_next = ll_dec;
                else
                    sol_f_next = interp1(pred_sorted,f_sorted,target_val(test_row),'linear','extrap'); % Determine inlet oxygen fraction for next segment using GPE results
                    if sol_f_next > ul_dec
                        sol_f_next = ul_dec;
                    elseif sol_f_next < ll_dec
                        sol_f_next = ll_dec;
                    else
                    end
                end

                % Adjust time for current segment if necessary to keep segment end target and actual value difference within oxygen concentration threshold
                t_vals = [test_dur(test_row), 200]; % Run past current segment to fully capture remnant increase
                f_vals = [test_f(test_row), sol_f_next];
                toggle = 0; % Define function output to be for increasing oxygen during segment
                target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                if abs(target_diff(1)) > O2_threshold || target_diff(2) ~= 1 % Difference between segment end target and actual value is not within oxygen concencentration threshold or there is no inflection point
                    if target_diff(1) > 0 || target_diff(2) ~= 1 % Actual value is lower than the target or there is not inflection --> increase run duration
                        target_reached = 0; % Initialize to enter while loop
                        counter_inc = 0;
                        while target_reached == 0
                            if target_diff(1) < 5 % Smaller step size if difference is within defined value
                                int = 1;
                            else
                                int = 5;
                            end
                            counter_inc = counter_inc + 1;
                            t_vals = [test_dur(test_row) + counter_inc*int, 200];
                            target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                            if abs(target_diff(1)) < O2_threshold
                                target_reached = 1; % End while loop
                            else
                                target_reached = 0; % Continue while loop
                            end
                            time_actual = test_dur(test_row) + counter_inc*int;
                            if time_actual > ul_time + 10
                                break % Exit while loop if run duration exceeds time threshold by 10 s
                            else
                            end
                        end
                    else % Actual value is higher than the target
                        target_reached = 0; % Initialize to enter while loop
                        counter_dec = 0;
                        while target_reached == 0
                            if target_diff(1) < 5 || time_actual < 5 % Smaller step size if difference is within defined value or adjusted time is less than 5 seconds
                                int = 1;
                            else
                                int = 5;
                            end
                            counter_dec = counter_dec + 1;
                            t_vals = [test_dur(test_row) - counter_dec*int, 200];
                            target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                            if abs(target_diff(1)) < O2_threshold
                                target_reached = 1; % End while loop
                            else
                                target_reached = 0; % Continue while loop
                            end
                            time_actual = test_dur(test_row) - counter_dec*int;
                            if time_actual < ll_time
                                break;
                            else
                            end
                        end
                    end
                elseif abs(target_diff(1)) < O2_threshold && target_diff(2) == 1 % Difference is within oxygen concentration threshold and inflection point is detected
                    time_actual = test_dur(test_row);
                else
                end
                sol_row = 1; % Solution identified; break out of while loop
            else % No solution found from GPE range comparison; continue with while loop
                sol_row = [];
                counter_loop = counter_loop + 1;
                if counter_loop > length(time_psol_sort)
                    break % No solution found; exit while loop
                else
                end
            end
        end

        % No solution found from GPE comparison
        if counter_loop > length(time_psol_sort)
            clear lms
            for n=1:length(time_psol_sort) % cycle through all possible solutions in filtered set
                test_row = find(time_val == time_psol_sort(n));
                test_row = test_row(1);
                test_x_1 = [ic.c3_20_ic*1e9, test_f(test_row), test_dur(test_row), ll_dec];
                pred_1 = predict(model_inc,test_x_1); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                test_x_2 = [ic.c3_20_ic*1e9, test_f(test_row), test_dur(test_row), ul_dec];
                pred_2 = predict(model_inc,test_x_2); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                lms(n) = (target_val(test_row) - pred_1)^2 + (target_val(test_row) - pred_2)^2;
            end
            min_lms = find(lms == min(lms)); % Find solution set with lowest difference between required and predicted remnant increase
            test_row = find(time_val == time_psol_sort(min_lms));
            test_row = test_row(1);
            test_x_1 = [ic.c3_20_ic*1e9, test_f(test_row), test_dur(test_row), ll_dec];
            pred_1 = predict(model_inc,test_x_1);  % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
            test_x_2 = [ic.c3_20_ic*1e9, test_f(test_row), test_dur(test_row), ul_dec];
            pred_2 = predict(model_inc,test_x_2);  % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
            test_x_3 = [ic.c3_20_ic*1e9, test_f(test_row), test_dur(test_row), mid_dec];
            pred_3 = predict(model_inc,test_x_3);  % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
            pred_range = [pred_1, pred_3, pred_2]; % GPE prediction range
            [pred_sorted, idx] = sort(pred_range);
            test_vals = [ll_dec, mid_dec, ul_dec];
            f_sorted = test_vals(idx);
            if length(unique(round(pred_range,3))) < length(pred_range) || range(pred_range) < range_threshold % There does not seem to be significant influence of inlet oxygen fraction for next segment
                sol_f_next = ll_dec;
            else
                sol_f_next = interp1(pred_sorted,f_sorted,target_val(test_row),'linear','extrap'); % Determine inlet oxygen fraction for next segment using GPE results
                if sol_f_next > ul_dec
                    sol_f_next = ul_dec;
                elseif sol_f_next < ll_dec
                    sol_f_next = ll_dec;
                else
                end
            end

            % Adjust time for current segment if necessary to keep segment end target and actual value difference within oxygen concentration threshold
            t_vals = [test_dur(test_row), 200]; % Run past current segment to fully capture remnant increase
            f_vals = [test_f(test_row), sol_f_next];
            toggle = 0; % Define function output to be for increasing segment
            target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
            if abs(target_diff(1)) > O2_threshold || target_diff(2) ~= 1 % Difference between segment end target and actual value is not within oxygen concencentration threshold or there is no inflection point
                if target_diff(1) > 0 || target_diff(2) ~= 1 % Actual value is lower than the target or there is not inflection --> increase run duration
                    target_reached = 0; % Initialize to enter while loop
                    counter_inc = 0;
                    while target_reached == 0
                        if target_diff(1) < 5 % Smaller step size if difference is within defined value
                            int = 1;
                        else
                            int = 5;
                        end
                        counter_inc = counter_inc + 1;
                        t_vals = [test_dur(test_row) + counter_inc*int, 200];
                        target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                        if abs(target_diff(1)) < O2_threshold
                            target_reached = 1; % End while loop
                        else
                            target_reached = 0; % Continue while loop
                        end
                        time_actual = test_dur(test_row) + counter_inc*int;
                        if time_actual > ul_time + 10
                            break;
                        else
                        end
                    end
                else % Actual value is higher than target
                    target_reached = 0; % Initialize to enter while loop
                    counter_dec = 0;
                    while target_reached == 0
                        if target_diff(1) < 5 || time_actual < 5 % Smaller step size if difference is within defined value or adjusted time is less than 5 seconds
                            int = 1;
                        else
                            int = 5;
                        end
                        counter_dec = counter_dec + 1;
                        t_vals = [test_dur(test_row) - counter_dec*int, 200];
                        target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                        if abs(target_diff(1)) < O2_threshold
                            target_reached = 1; % End while loop
                        else
                            target_reached = 0; % Continue while loop
                        end
                        time_actual = test_dur(test_row) - counter_dec*int;
                        if time_actual < ll_time
                            break;
                        else
                        end
                    end
                end
            elseif abs(target_diff(1)) < O2_threshold && target_diff(2) == 1 % Difference is within oxygen concentration threshold and inflection point is detected
                time_actual = test_dur(test_row);
            else
            end
        else % Solution already found
        end

        sol_f(i) = test_f(test_row);
        sol_dur(i) = time_actual;

    % Iteration > 2
    else
        sol_f(i) = sol_f_next; % Determined from previous iteration

        % Optimization targets for segment
        target_t = t1(indices(i+1)) - t1(indices(i)); % Target time for segment
        target_prev = target_c(i-1); % Target oxygen concentration at beginning of segment
        target = target_c(i); % Target oxygen concentration at end of segment

        % Run model with results of previous iteration to obtain the initial condition for current iteration
        f0 = sol_f(i-1); % Inlet oxygen fraction from previous iteration
        tdur = sol_dur(i-1); % Run duration from previous iteration
        toggle = 2; % Set function output to be the chamber state for initial condition
        sol_prev = opt_function_secondary(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,f0,tdur,target_prev,target,target_t,toggle,ic,O2_dir,T_in);
        ic = sol_prev; % Initial condition of chamber stored in structure

        % Primary optimization method
        clear objfun options
        % Determine whether segment will have increasing or decreasing oxygen concentration
        if target_c(i) > target_c(i-1) % Increasing oxygen concentration
            O2_dir = 0;
        elseif target_c(i) < target_c(i-1) % Decreasing oxygen concentration
            O2_dir = 1;
        else
        end

        % Adjust thresholds based on direction of oxygen profile for segment (time threshold is the same for both)
        if O2_dir == 0 % Increasing oxygen concentration
            lb = [ll_time, ll_dec];
            ub = [ul_time, ul_dec];
        else % Decreasing oxygen concentration
            lb = [ll_time, ll_inc];
            ub = [ul_time, ul_inc];
        end

        % Create optimization function and run Pareto front solver (gamultiobj)
        toggle = 1; % Define function output to be objective functions (toggle = 2 is to have function output be the chamber state for initial condition)
        objfun = @(fdur)opt_function_primary(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,sol_f(i),fdur(1),target_prev,target,target_t,toggle,ic,O2_dir,fdur(2),T_in);
        options = optimoptions('gamultiobj', 'Display','iter','MaxGenerations', 10,'PopulationSize',100,'ParetoFraction',0.5,'UseParallel', true);
        [sol, fval] = gamultiobj(objfun, 2, [], [], [], [], lb, ub, [], options); % Run optimization
        test_dur = sol(:,1); % Solutions for run duration of current segment
        test_next = sol(:,2); % Solutions for inlet oxygen fraction of previous segment

        % Results of objective functions
        time_val = fval(:,1);
        target_val = fval(:,2);
        target_prev_val = fval(:,3);

        % Find solution set from optimization results
        psol_ind = find(target_val < O2_threshold); % Filter possible solution set based on defined oxygen concentration threshold
        if double(isempty(psol_ind)) == 1
            psol_ind = 1; % Initialize to avoid code error
        else
        end
        time_psol = time_val(psol_ind);
        time_psol_sort = unique(sort(time_psol)); % Eliminate repeating solution sets
        sol_row = find(time_val == time_psol_sort(1)); % Select filtered solution set with smallest time difference
        if double(isempty(sol_row)) == 0
            sol_row = sol_row(1);
        else
        end

        if double(isempty(sol_row)) == 1 || time_val(sol_row) > time_threshold || target_val(sol_row) > O2_threshold % No solution found or solution set exceeds threshold limits
            disp('Primary method failed, running backup method...')

            % Backup method
            clear objfun options
            toggle = 1; % Define function output to be objective functions
            objfun = @(fdur)opt_function_secondary(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,sol_f(i),fdur,target_prev,target,target_t,toggle,ic,O2_dir,T_in);
            options = optimoptions('gamultiobj', 'Display','iter','MaxGenerations', 10,'PopulationSize',100,'ParetoFraction',0.5,'UseParallel', true);
            [sol, fval] = gamultiobj(objfun, 1, [], [], [], [], lb(1), ub(1), [], options); % Run optimization
            test_dur = sol; % Solutions for run duration (inlet oxygen fraction for current segment determined from previous iteration)

            % Results of objective functions
            time_val = fval(:,1);
            target_val = fval(:,2);
            target_prev_val = fval(:,3);

            % Find solution set from optimization results
            psol_ind = find(target_prev_val < O2_threshold); % Filter possible solutions based on defined oxygen concentration threshold
            time_psol = time_val(psol_ind);
            time_psol_sort = unique(sort(time_psol)); % Eliminate repeating solutions
            sol_row = []; % Initialize to enter while loop
            counter_loop = 1;
            if O2_dir == 0 % Increasing oxygen segment
                while double(isempty(sol_row)) == 1 % Run through possible solution sets until defined criteria is met
                    test_row = find(time_val == time_psol_sort(counter_loop));
                    test_row = test_row(1);
                    test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_dec];
                    pred_1 = predict(model_inc,test_x_1); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                    test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_dec];
                    pred_2 = predict(model_inc,test_x_2); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment

                    if target_val(test_row) > pred_1 - O2_threshold && target_val(test_row) < pred_2 + O2_threshold % Target remnant increase falls within GPE predicted range with uncertainty
                        test_x_3 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), mid_dec];
                        pred_3 = predict(model_inc,test_x_3); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                        pred_range = [pred_1, pred_3, pred_2]; % GPE prediction range
                        [pred_sorted, idx] = sort(pred_range);
                        test_vals = [ll_dec, mid_dec, ul_dec];
                        f_sorted = test_vals(idx);
                        if length(unique(round(pred_range,3))) < length(pred_range) || range(pred_range) < range_threshold % There does not seem to be significant influence of inlet oxygen fraction for next segment
                            sol_f_next = ll_dec;
                        else
                            sol_f_next = interp1(pred_sorted,f_sorted,target_val(test_row),'linear','extrap'); % Determine inlet oxygen fraction for next segment using GPE results
                            if sol_f_next > ul_dec
                                sol_f_next = ul_dec;
                            elseif sol_f_next < ll_dec
                                sol_f_next = ll_dec;
                            else
                            end
                        end

                        % Adjust time for current segment if necessary to keep segment end target and actual value difference within oxygen concentration threshold
                        t_vals = [test_dur(test_row), 200]; % Run past current segment to fully capture remnant increase
                        f_vals = [sol_f(i), sol_f_next];
                        toggle = 0; % Define function output to be for increasing oxygen during segment
                        target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                        if abs(target_diff(1)) > O2_threshold || target_diff(2) ~= 1 % Difference between segment end target and actual value is not within oxygen concentration threshold or there is no inflection point
                            if target_diff(1) > 0 || target_diff(2) ~= 1 % Actual value is lower than the target or there is not inflection --> increase run duration
                                target_reached = 0; % Initialize to enter while loop
                                counter_inc = 0;
                                while target_reached == 0
                                    if target_diff(1) < 5 % Smaller step size if difference is within defined value
                                        int = 1;
                                    else
                                        int = 5;
                                    end
                                    counter_inc = counter_inc + 1;
                                    t_vals = [test_dur(test_row) + counter_inc*int, 200];
                                    target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                    if abs(target_diff(1)) < O2_threshold 
                                        target_reached = 1; % End while loop
                                    else
                                        target_reached = 0; % Continue while loop
                                    end
                                    time_actual = test_dur(test_row) + counter_inc*int;
                                    if time_actual > ul_time + 10
                                        break % Exit while loop if run duration exceeds time threshold by 10 s
                                    else
                                    end
                                end
                            else % Actual value is higher than the target
                                target_reached = 0; % Initialize to enter while loop
                                counter_dec = 0;
                                while target_reached == 0
                                    if target_diff(1) < 5 || time_actual < 5 % Smaller step size if difference is within defined value or adjusted time is less than 5 seconds
                                        int = 1;
                                    else
                                        int = 5;
                                    end
                                    counter_dec = counter_dec + 1;
                                    t_vals = [test_dur(test_row) - counter_dec*int, 200];
                                    target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                    if abs(target_diff(1)) < O2_threshold
                                        target_reached = 1; % End while loop
                                    else
                                        target_reached = 0; % Continue while loop
                                    end
                                    time_actual = test_dur(test_row) - counter_dec*int;
                                    if time_actual < ll_time
                                        break;
                                    else
                                    end
                                end
                            end
                        elseif abs(target_diff(1)) < O2_threshold && target_diff(2) == 1 % Difference is within oxygen concentration threshold and inflection point is detected
                            time_actual = test_dur(test_row);
                        else
                        end
                        sol_row = 1; % Solution identified; break out of while loop
                    else % No solution found from GPE range comparison; continue with while loop
                        sol_row = [];
                        counter_loop = counter_loop + 1;
                        if counter_loop > length(time_psol_sort)
                            break % No solution found; exit while loop
                        else
                        end
                    end
                end

                % No solution found from GPE comparison
                if counter_loop > length(time_psol_sort)
                    clear lms
                    for n=1:length(time_psol_sort) % cycle through all possible solutions in filtered set
                        test_row = find(time_val == time_psol_sort(n));
                        test_row = test_row(1);
                        test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_dec];
                        pred_1 = predict(model_inc,test_x_1); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                        test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_dec];
                        pred_2 = predict(model_inc,test_x_2); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                        lms(n) = (target_val(test_row) - pred_1)^2 + (target_val(test_row) - pred_2)^2;
                    end
                    min_lms = find(lms == min(lms)); % Find solution set with lowest difference between required and predicted remnant increase
                    test_row = find(time_val == time_psol_sort(min_lms));
                    test_row = test_row(1);
                    test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_dec];
                    pred_1 = predict(model_inc,test_x_1);  % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                    test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_dec];
                    pred_2 = predict(model_inc,test_x_2);  % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                    test_x_3 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), mid_dec];
                    pred_3 = predict(model_inc,test_x_3);  % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                    pred_range = [pred_1, pred_3, pred_2]; % GPE prediction range
                    [pred_sorted, idx] = sort(pred_range);
                    test_vals = [ll_dec, mid_dec, ul_dec];
                    f_sorted = test_vals(idx);
                    if length(unique(round(pred_range,3))) < length(pred_range) || range(pred_range) < range_threshold % There does not seem to be significant influence of inlet oxygen fraction for next segment
                        sol_f_next = ll_dec;
                    else
                        sol_f_next = interp1(pred_sorted,f_sorted,target_val(test_row),'linear','extrap'); % Determine inlet oxygen fraction for next segment using GPE results
                        if sol_f_next > ul_dec
                            sol_f_next = ul_dec;
                        elseif sol_f_next < ll_dec
                            sol_f_next = ll_dec;
                        else
                        end
                    end

                    % Adjust time for current segment if necessary to keep segment end target and actual value difference within oxygen concentration threshold
                    t_vals = [test_dur(test_row), 200]; % Run past current segment to fully capture remnant increase
                    f_vals = [sol_f(i), sol_f_next];
                    toggle = 0; % Define function output to be for increasing segment
                    target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                    if abs(target_diff(1)) > O2_threshold || target_diff(2) ~= 1 % Difference between segment end target and actual value is not within oxygen concencentration threshold or there is no inflection point
                        if target_diff(1) > 0 || target_diff(2) ~= 1 % Actual value is lower than the target or there is not inflection --> increase run duration
                            target_reached = 0; % Initialize to enter while loop
                            counter_inc = 0;
                            while target_reached == 0
                                if target_diff(1) < 5 % Smaller step size if difference is within defined value
                                    int = 1;
                                else
                                    int = 5;
                                end
                                counter_inc = counter_inc + 1;
                                t_vals = [test_dur(test_row) + counter_inc*int, 200];
                                target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                if abs(target_diff(1)) < O2_threshold
                                    target_reached = 1; % End while loop
                                else
                                    target_reached = 0; % Continue while loop
                                end
                                time_actual = test_dur(test_row) + counter_inc*int;
                                if time_actual > ul_time + 10
                                    break % Exit while loop if run duration exceeds time threshold by 10 s
                                else
                                end
                            end
                        else % Actual value is higher than target
                            target_reached = 0; % Initialize to enter while loop
                            counter_dec = 0;
                            while target_reached == 0
                                if target_diff(1) < 5 || time_actual < 5 % Smaller step size if difference is within defined value or adjusted time is less than 5 seconds
                                    int = 1;
                                else
                                    int = 5;
                                end
                                counter_dec = counter_dec + 1;
                                t_vals = [test_dur(test_row) - counter_dec*int, 200];
                                target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                if abs(target_diff(1)) < O2_threshold
                                    target_reached = 1; % End while loop
                                else
                                    target_reached = 0; % Continue while loop
                                end
                                time_actual = test_dur(test_row) - counter_dec*int;
                                if time_actual < ll_time
                                    break;
                                else
                                end
                            end
                        end
                    elseif abs(target_diff(1)) < O2_threshold && target_diff(2) == 1 % Difference is within oxygen concentration threshold and inflection point is detected
                        time_actual = test_dur(test_row);
                    else
                    end
                else % Solution already found
                end
            else % Decreasing oxygen segment
                while double(isempty(sol_row)) == 1 % Run through possible solution sets until defined criteria is met
                    test_row = find(time_val == time_psol_sort(counter_loop));
                    test_row = test_row(1);
                    test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_inc];
                    pred_1 = predict(model_dec,test_x_1); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                    test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_inc];
                    pred_2 = predict(model_dec,test_x_2); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment

                    if target_val(test_row) > pred_2 - O2_threshold && target_val(test_row) < pred_1 + O2_threshold
                        test_x_3 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), mid_inc];
                        pred_3 = predict(model_dec,test_x_3); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                        pred_range = [pred_2, pred_3, pred_1]; % GPE prediction range
                        [pred_sorted, idx] = sort(pred_range);
                        test_vals = [ul_inc, mid_inc, ll_inc];
                        f_sorted = test_vals(idx);
                        if length(unique(round(pred_range,3))) < length(pred_range) || range(pred_range) < range_threshold % There does not seem to be significant influence of inlet oxygen fraction for next segment
                            sol_f_next = ul_inc;
                        else
                            sol_f_next = interp1(pred_sorted,f_sorted,target_val(test_row),'linear','extrap'); % Determine inlet oxygen fraction for next segment using GPE results
                            if sol_f_next > ul_inc
                                sol_f_next = ul_inc;
                            elseif sol_f_next < ll_inc
                                sol_f_next = ll_inc;
                            else
                            end
                        end

                        % Adjust time for current segment if necessary to keep segment end target and actual value difference within oxygen concentration threshold
                        t_vals = [test_dur(test_row), 200]; % Run past current segment to fully capture remnant decrease
                        f_vals = [sol_f(i), sol_f_next];
                        toggle = 1; % Define function output to be for decreasing oxygen during segment
                        target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                        if abs(target_diff(1)) > O2_threshold || target_diff(2) ~= 1  % Difference between segment end target and actual value is not within oxygen concencentration threshold or there is no inflection point
                            if target_diff(1) > 0 || target_diff(2) ~= 1 % Actual value is higher than target or there is no inflection --> increase run duration
                                target_reached = 0; % Initialize to enter while loop
                                counter_inc = 0;
                                while target_reached == 0 % Smaller step size if difference is within defined value
                                    if target_diff(1) < 5
                                        int = 1;
                                    else
                                        int = 5;
                                    end
                                    counter_inc = counter_inc + 1;
                                    t_vals = [test_dur(test_row) + counter_inc*int, 200];
                                    target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                    if abs(target_diff(1)) < O2_threshold
                                        target_reached = 1; % End while loop
                                    else
                                        target_reached = 0; % Continue while loop
                                    end
                                    time_actual = test_dur(test_row) + counter_inc*int;
                                    if time_actual > ul_time + 10
                                        break % Exit while loop if run duration exceeds time threshold by 10 s
                                    else
                                    end
                                end
                            else % Actual value is lower than target
                                target_reached = 0; % Initialize to enter while loop
                                counter_dec = 0;
                                while target_reached == 0
                                    if target_diff(1) < 5 || time_actual < 5 % Smaller step size if difference is within defined value or adjusted time is less than 5 seconds
                                        int = 1;
                                    else
                                        int = 5;
                                    end
                                    counter_dec = counter_dec + 1;
                                    t_vals = [test_dur(test_row) - counter_dec*int, 200];
                                    target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                    if abs(target_diff(1)) < O2_threshold
                                        target_reached = 1; % End while loop
                                    else
                                        target_reached = 0; % Continue while loop
                                    end
                                    time_actual = test_dur(test_row) - counter_dec*int;
                                    if time_actual < ll_time
                                        break;
                                    else
                                    end
                                end
                            end
                        elseif abs(target_diff(1)) < O2_threshold && target_diff(2) == 1 % Difference is within oxygen concentration threshold and inflection point is detected
                            time_actual = test_dur(test_row);
                        else
                        end
                        sol_row = 1; % Solution identified; break out of while loop
                    else % No solution found from GPE range comparison; continue with while loop
                        sol_row = [];
                        counter_loop = counter_loop + 1;
                        if counter_loop > length(time_psol_sort)
                            break % No solution found; exit while loop
                        else
                        end
                    end
                end

                % No solution found from GPE comparison
                if counter_loop > length(time_psol_sort) % Cycle through all possible solutions in filtered set
                    clear lms
                    for n=1:length(time_psol_sort)
                        test_row = find(time_val == time_psol_sort(n));
                        test_row = test_row(1);
                        test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_inc];
                        pred_1 = predict(model_inc,test_x_1); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                        test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_inc];
                        pred_2 = predict(model_inc,test_x_2); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                        lms(n) = (target_val(test_row) - pred_1)^2 + (target_val(test_row) - pred_2)^2;
                    end
                    min_lms = find(lms == min(lms)); % Find solution set with lowest difference between required and predicted remnant increase
                    test_row = find(time_val == time_psol_sort(min_lms));
                    test_row = test_row(1);
                    test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_inc];
                    pred_1 = predict(model_inc,test_x_1); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                    test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_inc];
                    pred_2 = predict(model_inc,test_x_2); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                    test_x_3 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), mid_inc];
                    pred_3 = predict(model_inc,test_x_3); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                    pred_range = [pred_2, pred_3, pred_1]; % GPE prediction range
                    [pred_sorted, idx] = sort(pred_range);
                    test_vals = [ul_inc, mid_inc, ll_inc];
                    f_sorted = test_vals(idx);
                    if length(unique(round(pred_range,3))) < length(pred_range) || range(pred_range) < range_threshold % There does not seem to be significant influence of inlet oxygen fraction for next segment
                        sol_f_next = ul_inc;
                    else
                        sol_f_next = interp1(pred_sorted,f_sorted,target_val(test_row),'linear','extrap'); % Determine inlet oxygen fraction for next segment using GPE results
                        if sol_f_next > ul_inc
                            sol_f_next = ul_inc;
                        elseif sol_f_next < ll_inc
                            sol_f_next = ll_inc;
                        else
                        end
                    end

                    % Adjust time for current segment if necessary to keep segment end target and actual value difference within oxygen concentration threshold
                    t_vals = [test_dur(test_row), 200]; % Run past current segment to fully capture remnant increase
                    f_vals = [sol_f(i), sol_f_next];
                    toggle = 1; % Define function output to be for decreasing segment
                    target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                    if abs(target_diff(1)) > O2_threshold || target_diff(2) ~= 1 % Difference between segment end target and actual value is not within oxygen concentration threshold or there is no inflection point
                        if target_diff(1) > 0 || target_diff(2) ~= 1
                            target_reached = 0; % Initialize to enter while loop
                            counter_inc = 0;
                            while target_reached == 0
                                if target_diff(1) < 5 % Smaller step size if difference is within defined value
                                    int = 1;
                                else
                                    int = 5;
                                end
                                counter_inc = counter_inc + 1;
                                t_vals = [test_dur(test_row) + counter_inc*int, 200];
                                target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                if abs(target_diff(1)) < O2_threshold
                                    target_reached = 1; % End while loop
                                else
                                    target_reached = 0; % Continue while loop
                                end
                                time_actual = test_dur(test_row) + counter_inc*int;
                                if time_actual > ul_time + 10
                                    break % Exit while loop if run duration exceeds time threshold by 10 s
                                else
                                end
                            end
                        else % Actual value is lower than target
                            target_reached = 0; % Initialize to enter while loop
                            counter_dec = 0;
                            while target_reached == 0
                                if target_diff(1) < 5 || time_actual < 5 % Smaller step size if difference is within defined value or adjusted time is less than 5 seconds
                                    int = 1;
                                else
                                    int = 5;
                                end
                                counter_dec = counter_dec + 1;
                                t_vals = [test_dur(test_row) - counter_dec*int, 200];
                                target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                if abs(target_diff(1)) < O2_threshold
                                    target_reached = 1; % End while loop
                                else
                                    target_reached = 0; % Continue while loop
                                end
                                time_actual = test_dur(test_row) - counter_dec*int;
                                if time_actual < ll_time
                                    break;
                                else
                                end
                            end
                        end
                    elseif abs(target_diff(1)) < O2_threshold && target_diff(2) == 1 % Difference is within oxygen concentration threshold and inflection point is detected
                        time_actual = test_dur(test_row);
                    else
                    end
                else % Solution already found
                end
            end
            sol_dur(i) = time_actual;

            % Run duration exceeds upper limit or is too short (low amplitude for target segment)
            if target_diff(3) > time_threshold
                time_dir = 1;
            elseif target_diff(3) < -time_threshold
                time_dir = 2;
            elseif target_diff(3) == 100
                time_dir = 3;
            else
                time_dir = 0;
            end

            % If time of segment falls outside threshold, adjust inlet oxygen fraction and re-run backup method
            if time_dir >= 1
                if time_dir == 3
                    if O2_dir == 0 % For increasing oxygen segment
                        sol_f(i) = ul_inc;
                    else % For decreasing oxygen segment
                        sol_f(i) = ll_dec;
                    end
                elseif time_dir == 1
                    % Re-run backup method
                    while target_diff(3) > time_threshold
                        if O2_dir == 0
                            sol_f(i) = sol_f(i) - O2f_step;
                        else
                            sol_f(i) = sol_f(i) + O2f_step;
                        end

                        if O2_dir == 0
                            if sol_f(i) <= ll_inc
                                sol_f(i) = sol_f(i) + O2f_step;
                                break
                            end
                        else
                            if sol_f(i) >= ul_dec
                                sol_f(i) = sol_f(i) - O2f_step;
                                break
                            end
                        end

                        clear objfun options
                        toggle = 1; % Define function output to be objective functions
                        objfun = @(fdur)opt_function_secondary(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,sol_f(i),fdur,target_prev,target,target_t,toggle,ic,O2_dir,T_in);
                        options = optimoptions('gamultiobj', 'Display','iter','MaxGenerations', 10,'PopulationSize',100,'ParetoFraction',0.5,'UseParallel', true);
                        [sol, fval] = gamultiobj(objfun, 1, [], [], [], [], lb(1), ub(1), [], options); % Run optimization
                        test_dur = sol; % Solutions for run duration (inlet oxygen fraction for current segment determined from previous iteration)

                        % Results of objective functions
                        time_val = fval(:,1);
                        target_val = fval(:,2);
                        target_prev_val = fval(:,3);

                        % Find solution set from optimization results
                        psol_ind = find(target_prev_val < O2_threshold); % Filter possible solutions based on defined oxygen concentration threshold
                        time_psol = time_val(psol_ind);
                        time_psol_sort = unique(sort(time_psol)); % Eliminate repeating solutions
                        sol_row = []; % Initialize to enter while loop
                        counter_loop = 1;
                        if O2_dir == 0 % Increasing oxygen segment
                            while double(isempty(sol_row)) == 1 % Run through possible solution sets until defined criteria is met
                                test_row = find(time_val == time_psol_sort(counter_loop));
                                test_row = test_row(1);
                                test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_dec];
                                pred_1 = predict(model_inc,test_x_1); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_dec];
                                pred_2 = predict(model_inc,test_x_2); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment

                                if target_val(test_row) > pred_1 - O2_threshold && target_val(test_row) < pred_2 + O2_threshold % Target remnant increase falls within GPE predicted range with uncertainty
                                    test_x_3 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), mid_dec];
                                    pred_3 = predict(model_inc,test_x_3); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                    pred_range = [pred_1, pred_3, pred_2]; % GPE prediction range
                                    [pred_sorted, idx] = sort(pred_range);
                                    test_vals = [ll_dec, mid_dec, ul_dec];
                                    f_sorted = test_vals(idx);
                                    if length(unique(round(pred_range,3))) < length(pred_range) || range(pred_range) < range_threshold % There does not seem to be significant influence of inlet oxygen fraction for next segment
                                        sol_f_next = ll_dec;
                                    else
                                        sol_f_next = interp1(pred_sorted,f_sorted,target_val(test_row),'linear','extrap'); % Determine inlet oxygen fraction for next segment using GPE results
                                        if sol_f_next > ul_dec
                                            sol_f_next = ul_dec;
                                        elseif sol_f_next < ll_dec
                                            sol_f_next = ll_dec;
                                        else
                                        end
                                    end

                                    % Adjust time for current segment if necessary to keep segment end target and actual value difference within oxygen concentration threshold
                                    t_vals = [test_dur(test_row), 200]; % Run past current segment to fully capture remnant increase
                                    f_vals = [sol_f(i), sol_f_next];
                                    toggle = 0; % Define function output to be for increasing oxygen during segment
                                    target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                    if abs(target_diff(1)) > O2_threshold || target_diff(2) ~= 1 % Difference between segment end target and actual value is not within oxygen concencentration threshold or there is no inflection point
                                        if target_diff(1) > 0 || target_diff(2) ~= 1 % Actual value is lower than the target or there is not inflection --> increase run duration
                                            target_reached = 0; % Initialize to enter while loop
                                            counter_inc = 0;
                                            while target_reached == 0
                                                if target_diff(1) < 5 % Smaller step size if difference is within defined value
                                                    int = 1;
                                                else
                                                    int = 5;
                                                end
                                                counter_inc = counter_inc + 1;
                                                t_vals = [test_dur(test_row) + counter_inc*int, 200];
                                                target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                                if abs(target_diff(1)) < O2_threshold
                                                    target_reached = 1; % End while loop
                                                else
                                                    target_reached = 0; % Continue while loop
                                                end
                                                time_actual = test_dur(test_row) + counter_inc*int;
                                                if time_actual > ul_time + 10
                                                    break % Exit while loop if run duration exceeds time threshold by 10 s
                                                else
                                                end
                                            end
                                        else % Actual value is higher than the target
                                            target_reached = 0; % Initialize to enter while loop
                                            counter_dec = 0;
                                            while target_reached == 0
                                                if target_diff(1) < 5 || time_actual < 5 % Smaller step size if difference is within defined value or adjusted time is less than 5 seconds
                                                    int = 1;
                                                else
                                                    int = 5;
                                                end
                                                counter_dec = counter_dec + 1;
                                                t_vals = [test_dur(test_row) - counter_dec*int, 200];
                                                target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                                if abs(target_diff(1)) < O2_threshold
                                                    target_reached = 1; % End while loop
                                                else
                                                    target_reached = 0; % Continue while loop
                                                end
                                                time_actual = test_dur(test_row) - counter_dec*int;
                                                if time_actual < ll_time
                                                    break;
                                                else
                                                end
                                            end
                                        end
                                    elseif abs(target_diff(1)) < O2_threshold && target_diff(2) == 1 % Difference is within oxygen concentration threshold and inflection point is detected
                                        time_actual = test_dur(test_row);
                                    else
                                    end
                                    sol_row = 1; % Solution identified; break out of while loop
                                else % No solution found from GPE range comparison; continue with while loop
                                    sol_row = [];
                                    counter_loop = counter_loop + 1;
                                    if counter_loop > length(time_psol_sort)
                                        break % No solution found; exit while loop
                                    else
                                    end
                                end
                            end

                            % No solution found from GPE comparison
                            if counter_loop > length(time_psol_sort)
                                clear lms
                                for n=1:length(time_psol_sort) % cycle through all possible solutions in filtered set
                                    test_row = find(time_val == time_psol_sort(n));
                                    test_row = test_row(1);
                                    test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_dec];
                                    pred_1 = predict(model_inc,test_x_1); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                    test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_dec];
                                    pred_2 = predict(model_inc,test_x_2); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                    lms(n) = (target_val(test_row) - pred_1)^2 + (target_val(test_row) - pred_2)^2;
                                end
                                min_lms = find(lms == min(lms)); % Find solution set with lowest difference between required and predicted remnant increase
                                test_row = find(time_val == time_psol_sort(min_lms));
                                test_row = test_row(1);
                                test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_dec];
                                pred_1 = predict(model_inc,test_x_1);  % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_dec];
                                pred_2 = predict(model_inc,test_x_2);  % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                test_x_3 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), mid_dec];
                                pred_3 = predict(model_inc,test_x_3);  % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                pred_range = [pred_1, pred_3, pred_2]; % GPE prediction range
                                [pred_sorted, idx] = sort(pred_range);
                                test_vals = [ll_dec, mid_dec, ul_dec];
                                f_sorted = test_vals(idx);
                                if length(unique(round(pred_range,3))) < length(pred_range) || range(pred_range) < range_threshold % There does not seem to be significant influence of inlet oxygen fraction for next segment
                                    sol_f_next = ll_dec;
                                else
                                    sol_f_next = interp1(pred_sorted,f_sorted,target_val(test_row),'linear','extrap'); % Determine inlet oxygen fraction for next segment using GPE results
                                    if sol_f_next > ul_dec
                                        sol_f_next = ul_dec;
                                    elseif sol_f_next < ll_dec
                                        sol_f_next = ll_dec;
                                    else
                                    end
                                end

                                % Adjust time for current segment if necessary to keep segment end target and actual value difference within oxygen concentration threshold
                                t_vals = [test_dur(test_row), 200]; % Run past current segment to fully capture remnant increase
                                f_vals = [sol_f(i), sol_f_next];
                                toggle = 0; % Define function output to be for increasing segment
                                target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                if abs(target_diff(1)) > O2_threshold || target_diff(2) ~= 1 % Difference between segment end target and actual value is not within oxygen concencentration threshold or there is no inflection point
                                    if target_diff(1) > 0 || target_diff(2) ~= 1 % Actual value is lower than the target or there is not inflection --> increase run duration
                                        target_reached = 0; % Initialize to enter while loop
                                        counter_inc = 0;
                                        while target_reached == 0
                                            if target_diff(1) < 5 % Smaller step size if difference is within defined value
                                                int = 1;
                                            else
                                                int = 5;
                                            end
                                            counter_inc = counter_inc + 1;
                                            t_vals = [test_dur(test_row) + counter_inc*int, 200];
                                            target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                            if abs(target_diff(1)) < O2_threshold
                                                target_reached = 1; % End while loop
                                            else
                                                target_reached = 0; % Continue while loop
                                            end
                                            time_actual = test_dur(test_row) + counter_inc*int;
                                            if time_actual > ul_time + 10
                                                break % Exit while loop if run duration exceeds time threshold by 10 s
                                            else
                                            end
                                        end
                                    else % Actual value is higher than target
                                        target_reached = 0; % Initialize to enter while loop
                                        counter_dec = 0;
                                        while target_reached == 0
                                            if target_diff(1) < 5 || time_actual < 5 % Smaller step size if difference is within defined value or adjusted time is less than 5 seconds
                                                int = 1;
                                            else
                                                int = 5;
                                            end
                                            counter_dec = counter_dec + 1;
                                            t_vals = [test_dur(test_row) - counter_dec*int, 200];
                                            target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                            if abs(target_diff(1)) < O2_threshold
                                                target_reached = 1; % End while loop
                                            else
                                                target_reached = 0; % Continue while loop
                                            end
                                            time_actual = test_dur(test_row) - counter_dec*int;
                                            if time_actual < ll_time
                                                break;
                                            else
                                            end
                                        end
                                    end
                                elseif abs(target_diff(1)) < O2_threshold && target_diff(2) == 1 % Difference is within oxygen concentration threshold and inflection point is detected
                                    time_actual = test_dur(test_row);
                                else
                                end
                            else % Solution already found
                            end
                        else % Decreasing oxygen segment
                            while double(isempty(sol_row)) == 1 % Run through possible solution sets until defined criteria is met
                                test_row = find(time_val == time_psol_sort(counter_loop));
                                test_row = test_row(1);
                                test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_inc];
                                pred_1 = predict(model_dec,test_x_1); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_inc];
                                pred_2 = predict(model_dec,test_x_2); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment

                                if target_val(test_row) > pred_2 - O2_threshold && target_val(test_row) < pred_1 + O2_threshold
                                    test_x_3 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), mid_inc];
                                    pred_3 = predict(model_dec,test_x_3); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                    pred_range = [pred_2, pred_3, pred_1]; % GPE prediction range
                                    [pred_sorted, idx] = sort(pred_range);
                                    test_vals = [ul_inc, mid_inc, ll_inc];
                                    f_sorted = test_vals(idx);
                                    if length(unique(round(pred_range,3))) < length(pred_range) || range(pred_range) < range_threshold % There does not seem to be significant influence of inlet oxygen fraction for next segment
                                        sol_f_next = ul_inc;
                                    else
                                        sol_f_next = interp1(pred_sorted,f_sorted,target_val(test_row),'linear','extrap'); % Determine inlet oxygen fraction for next segment using GPE results
                                        if sol_f_next > ul_inc
                                            sol_f_next = ul_inc;
                                        elseif sol_f_next < ll_inc
                                            sol_f_next = ll_inc;
                                        else
                                        end
                                    end

                                    % Adjust time for current segment if necessary to keep segment end target and actual value difference within oxygen concentration threshold
                                    t_vals = [test_dur(test_row), 200]; % Run past current segment to fully capture remnant decrease
                                    f_vals = [sol_f(i), sol_f_next];
                                    toggle = 1; % Define function output to be for decreasing oxygen during segment
                                    target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                    if abs(target_diff(1)) > O2_threshold || target_diff(2) ~= 1  % Difference between segment end target and actual value is not within oxygen concencentration threshold or there is no inflection point
                                        if target_diff(1) > 0 || target_diff(2) ~= 1 % Actual value is higher than target or there is no inflection --> increase run duration
                                            target_reached = 0; % Initialize to enter while loop
                                            counter_inc = 0;
                                            while target_reached == 0 % Smaller step size if difference is within defined value
                                                if target_diff(1) < 5
                                                    int = 1;
                                                else
                                                    int = 5;
                                                end
                                                counter_inc = counter_inc + 1;
                                                t_vals = [test_dur(test_row) + counter_inc*int, 200];
                                                target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                                if abs(target_diff(1)) < O2_threshold
                                                    target_reached = 1; % End while loop
                                                else
                                                    target_reached = 0; % Continue while loop
                                                end
                                                time_actual = test_dur(test_row) + counter_inc*int;
                                                if time_actual > ul_time + 10
                                                    break % Exit while loop if run duration exceeds time threshold by 10 s
                                                else
                                                end
                                            end
                                        else % Actual value is lower than target
                                            target_reached = 0; % Initialize to enter while loop
                                            counter_dec = 0;
                                            while target_reached == 0
                                                if target_diff(1) < 5 || time_actual < 5 % Smaller step size if difference is within defined value or adjusted time is less than 5 seconds
                                                    int = 1;
                                                else
                                                    int = 5;
                                                end
                                                counter_dec = counter_dec + 1;
                                                t_vals = [test_dur(test_row) - counter_dec*int, 200];
                                                target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                                if abs(target_diff(1)) < O2_threshold
                                                    target_reached = 1; % End while loop
                                                else
                                                    target_reached = 0; % Continue while loop
                                                end
                                                time_actual = test_dur(test_row) - counter_dec*int;
                                                if time_actual < ll_time
                                                    break;
                                                else
                                                end
                                            end
                                        end
                                    elseif abs(target_diff(1)) < O2_threshold && target_diff(2) == 1 % Difference is within oxygen concentration threshold and inflection point is detected
                                        time_actual = test_dur(test_row);
                                    else
                                    end
                                    sol_row = 1; % Solution identified; break out of while loop
                                else % No solution found from GPE range comparison; continue with while loop
                                    sol_row = [];
                                    counter_loop = counter_loop + 1;
                                    if counter_loop > length(time_psol_sort)
                                        break % No solution found; exit while loop
                                    else
                                    end
                                end
                            end

                            % No solution found from GPE comparison
                            if counter_loop > length(time_psol_sort) % Cycle through all possible solutions in filtered set
                                clear lms
                                for n=1:length(time_psol_sort)
                                    test_row = find(time_val == time_psol_sort(n));
                                    test_row = test_row(1);
                                    test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_inc];
                                    pred_1 = predict(model_inc,test_x_1); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                    test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_inc];
                                    pred_2 = predict(model_inc,test_x_2); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                    lms(n) = (target_val(test_row) - pred_1)^2 + (target_val(test_row) - pred_2)^2;
                                end
                                min_lms = find(lms == min(lms)); % Find solution set with lowest difference between required and predicted remnant increase
                                test_row = find(time_val == time_psol_sort(min_lms));
                                test_row = test_row(1);
                                test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_inc];
                                pred_1 = predict(model_inc,test_x_1); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_inc];
                                pred_2 = predict(model_inc,test_x_2); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                test_x_3 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), mid_inc];
                                pred_3 = predict(model_inc,test_x_3); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                pred_range = [pred_2, pred_3, pred_1]; % GPE prediction range
                                [pred_sorted, idx] = sort(pred_range);
                                test_vals = [ul_inc, mid_inc, ll_inc];
                                f_sorted = test_vals(idx);
                                if length(unique(round(pred_range,3))) < length(pred_range) || range(pred_range) < range_threshold % There does not seem to be significant influence of inlet oxygen fraction for next segment
                                    sol_f_next = ul_inc;
                                else
                                    sol_f_next = interp1(pred_sorted,f_sorted,target_val(test_row),'linear','extrap'); % Determine inlet oxygen fraction for next segment using GPE results
                                    if sol_f_next > ul_inc
                                        sol_f_next = ul_inc;
                                    elseif sol_f_next < ll_inc
                                        sol_f_next = ll_inc;
                                    else
                                    end
                                end

                                % Adjust time for current segment if necessary to keep segment end target and actual value difference within oxygen concentration threshold
                                t_vals = [test_dur(test_row), 200]; % Run past current segment to fully capture remnant increase
                                f_vals = [sol_f(i), sol_f_next];
                                toggle = 1; % Define function output to be for decreasing segment
                                target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                if abs(target_diff(1)) > O2_threshold || target_diff(2) ~= 1 % Difference between segment end target and actual value is not within oxygen concentration threshold or there is no inflection point
                                    if target_diff(1) > 0 || target_diff(2) ~= 1
                                        target_reached = 0; % Initialize to enter while loop
                                        counter_inc = 0;
                                        while target_reached == 0
                                            if target_diff(1) < 5 % Smaller step size if difference is within defined value
                                                int = 1;
                                            else
                                                int = 5;
                                            end
                                            counter_inc = counter_inc + 1;
                                            t_vals = [test_dur(test_row) + counter_inc*int, 200];
                                            target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                            if abs(target_diff(1)) < O2_threshold
                                                target_reached = 1; % End while loop
                                            else
                                                target_reached = 0; % Continue while loop
                                            end
                                            time_actual = test_dur(test_row) + counter_inc*int;
                                            if time_actual > ul_time + 10
                                                break % Exit while loop if run duration exceeds time threshold by 10 s
                                            else
                                            end
                                        end
                                    else % Actual value is lower than target
                                        target_reached = 0; % Initialize to enter while loop
                                        counter_dec = 0;
                                        while target_reached == 0
                                            if target_diff(1) < 5 || time_actual < 5 % Smaller step size if difference is within defined value or adjusted time is less than 5 seconds
                                                int = 1;
                                            else
                                                int = 5;
                                            end
                                            counter_dec = counter_dec + 1;
                                            t_vals = [test_dur(test_row) - counter_dec*int, 200];
                                            target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                            if abs(target_diff(1)) < O2_threshold
                                                target_reached = 1; % End while loop
                                            else
                                                target_reached = 0; % Continue while loop
                                            end
                                            time_actual = test_dur(test_row) - counter_dec*int;
                                            if time_actual < ll_time
                                                break;
                                            else
                                            end
                                        end
                                    end
                                elseif abs(target_diff(1)) < O2_threshold && target_diff(2) == 1 % Difference is within oxygen concentration threshold and inflection point is detected
                                    time_actual = test_dur(test_row);
                                else
                                end
                            else % Solution already found
                            end
                        end
                    end
                    sol_dur(i) = time_actual;
                elseif time_dir == 2
                    % Re-run backup method
                    while target_diff(3) < -time_threshold
                        if O2_dir == 0
                            sol_f(i) = sol_f(i) + O2f_step;
                        else
                            sol_f(i) = sol_f(i) - O2f_step;
                        end

                        if O2_dir == 0
                            if sol_f(i) >= ul_inc
                                sol_f(i) = sol_f(i) - O2f_step;
                                break
                            end
                        else
                            if sol_f(i) <= ll_dec
                                sol_f(i) = sol_f(i) + O2f_step;
                                break
                            end
                        end
                        
                        clear objfun options
                        toggle = 1; % Define function output to be objective functions
                        objfun = @(fdur)opt_function_secondary(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,sol_f(i),fdur,target_prev,target,target_t,toggle,ic,O2_dir,T_in);
                        options = optimoptions('gamultiobj', 'Display','iter','MaxGenerations', 10,'PopulationSize',100,'ParetoFraction',0.5,'UseParallel', true);
                        [sol, fval] = gamultiobj(objfun, 1, [], [], [], [], lb(1), ub(1), [], options); % Run optimization
                        test_dur = sol; % Solutions for run duration (inlet oxygen fraction for current segment determined from previous iteration)

                        % Results of objective functions
                        time_val = fval(:,1);
                        target_val = fval(:,2);
                        target_prev_val = fval(:,3);

                        % Find solution set from optimization results
                        psol_ind = find(target_prev_val < O2_threshold); % Filter possible solutions based on defined oxygen concentration threshold
                        time_psol = time_val(psol_ind);
                        time_psol_sort = unique(sort(time_psol)); % Eliminate repeating solutions
                        sol_row = []; % Initialize to enter while loop
                        counter_loop = 1;
                        if O2_dir == 0 % Increasing oxygen segment
                            while double(isempty(sol_row)) == 1 % Run through possible solution sets until defined criteria is met
                                test_row = find(time_val == time_psol_sort(counter_loop));
                                test_row = test_row(1);
                                test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_dec];
                                pred_1 = predict(model_inc,test_x_1); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_dec];
                                pred_2 = predict(model_inc,test_x_2); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment

                                if target_val(test_row) > pred_1 - O2_threshold && target_val(test_row) < pred_2 + O2_threshold % Target remnant increase falls within GPE predicted range with uncertainty
                                    test_x_3 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), mid_dec];
                                    pred_3 = predict(model_inc,test_x_3); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                    pred_range = [pred_1, pred_3, pred_2]; % GPE prediction range
                                    [pred_sorted, idx] = sort(pred_range);
                                    test_vals = [ll_dec, mid_dec, ul_dec];
                                    f_sorted = test_vals(idx);
                                    if length(unique(round(pred_range,3))) < length(pred_range) || range(pred_range) < range_threshold % There does not seem to be significant influence of inlet oxygen fraction for next segment
                                        sol_f_next = ll_dec;
                                    else
                                        sol_f_next = interp1(pred_sorted,f_sorted,target_val(test_row),'linear','extrap'); % Determine inlet oxygen fraction for next segment using GPE results
                                        if sol_f_next > ul_dec
                                            sol_f_next = ul_dec;
                                        elseif sol_f_next < ll_dec
                                            sol_f_next = ll_dec;
                                        else
                                        end
                                    end

                                    % Adjust time for current segment if necessary to keep segment end target and actual value difference within oxygen concentration threshold
                                    t_vals = [test_dur(test_row), 200]; % Run past current segment to fully capture remnant increase
                                    f_vals = [sol_f(i), sol_f_next];
                                    toggle = 0; % Define function output to be for increasing oxygen during segment
                                    target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                    if abs(target_diff(1)) > O2_threshold || target_diff(2) ~= 1 % Difference between segment end target and actual value is not within oxygen concencentration threshold or there is no inflection point
                                        if target_diff(1) > 0 || target_diff(2) ~= 1 % Actual value is lower than the target or there is not inflection --> increase run duration
                                            target_reached = 0; % Initialize to enter while loop
                                            counter_inc = 0;
                                            while target_reached == 0
                                                if target_diff(1) < 5 % Smaller step size if difference is within defined value
                                                    int = 1;
                                                else
                                                    int = 5;
                                                end
                                                counter_inc = counter_inc + 1;
                                                t_vals = [test_dur(test_row) + counter_inc*int, 200];
                                                target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                                if abs(target_diff(1)) < O2_threshold
                                                    target_reached = 1; % End while loop
                                                else
                                                    target_reached = 0; % Continue while loop
                                                end
                                                time_actual = test_dur(test_row) + counter_inc*int;
                                                if time_actual > ul_time + 10
                                                    break % Exit while loop if run duration exceeds time threshold by 10 s
                                                else
                                                end
                                            end
                                        else % Actual value is higher than the target
                                            target_reached = 0; % Initialize to enter while loop
                                            counter_dec = 0;
                                            while target_reached == 0
                                                if target_diff(1) < 5 || time_actual < 5 % Smaller step size if difference is within defined value or adjusted time is less than 5 seconds
                                                    int = 1;
                                                else
                                                    int = 5;
                                                end
                                                counter_dec = counter_dec + 1;
                                                t_vals = [test_dur(test_row) - counter_dec*int, 200];
                                                target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                                if abs(target_diff(1)) < O2_threshold
                                                    target_reached = 1; % End while loop
                                                else
                                                    target_reached = 0; % Continue while loop
                                                end
                                                time_actual = test_dur(test_row) - counter_dec*int;
                                                if time_actual < ll_time
                                                    break;
                                                else
                                                end
                                            end
                                        end
                                    elseif abs(target_diff(1)) < O2_threshold && target_diff(2) == 1 % Difference is within oxygen concentration threshold and inflection point is detected
                                        time_actual = test_dur(test_row);
                                    else
                                    end
                                    sol_row = 1; % Solution identified; break out of while loop
                                else % No solution found from GPE range comparison; continue with while loop
                                    sol_row = [];
                                    counter_loop = counter_loop + 1;
                                    if counter_loop > length(time_psol_sort)
                                        break % No solution found; exit while loop
                                    else
                                    end
                                end
                            end

                            % No solution found from GPE comparison
                            if counter_loop > length(time_psol_sort)
                                clear lms
                                for n=1:length(time_psol_sort) % cycle through all possible solutions in filtered set
                                    test_row = find(time_val == time_psol_sort(n));
                                    test_row = test_row(1);
                                    test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_dec];
                                    pred_1 = predict(model_inc,test_x_1); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                    test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_dec];
                                    pred_2 = predict(model_inc,test_x_2); % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                    lms(n) = (target_val(test_row) - pred_1)^2 + (target_val(test_row) - pred_2)^2;
                                end
                                min_lms = find(lms == min(lms)); % Find solution set with lowest difference between required and predicted remnant increase
                                test_row = find(time_val == time_psol_sort(min_lms));
                                test_row = test_row(1);
                                test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_dec];
                                pred_1 = predict(model_inc,test_x_1);  % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_dec];
                                pred_2 = predict(model_inc,test_x_2);  % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                test_x_3 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), mid_dec];
                                pred_3 = predict(model_inc,test_x_3);  % Predict remnant increase using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                pred_range = [pred_1, pred_3, pred_2]; % GPE prediction range
                                [pred_sorted, idx] = sort(pred_range);
                                test_vals = [ll_dec, mid_dec, ul_dec];
                                f_sorted = test_vals(idx);
                                if length(unique(round(pred_range,3))) < length(pred_range) || range(pred_range) < range_threshold % There does not seem to be significant influence of inlet oxygen fraction for next segment
                                    sol_f_next = ll_dec;
                                else
                                    sol_f_next = interp1(pred_sorted,f_sorted,target_val(test_row),'linear','extrap'); % Determine inlet oxygen fraction for next segment using GPE results
                                    if sol_f_next > ul_dec
                                        sol_f_next = ul_dec;
                                    elseif sol_f_next < ll_dec
                                        sol_f_next = ll_dec;
                                    else
                                    end
                                end

                                % Adjust time for current segment if necessary to keep segment end target and actual value difference within oxygen concentration threshold
                                t_vals = [test_dur(test_row), 200]; % Run past current segment to fully capture remnant increase
                                f_vals = [sol_f(i), sol_f_next];
                                toggle = 0; % Define function output to be for increasing segment
                                target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                if abs(target_diff(1)) > O2_threshold || target_diff(2) ~= 1 % Difference between segment end target and actual value is not within oxygen concencentration threshold or there is no inflection point
                                    if target_diff(1) > 0 || target_diff(2) ~= 1 % Actual value is lower than the target or there is not inflection --> increase run duration
                                        target_reached = 0; % Initialize to enter while loop
                                        counter_inc = 0;
                                        while target_reached == 0
                                            if target_diff(1) < 5 % Smaller step size if difference is within defined value
                                                int = 1;
                                            else
                                                int = 5;
                                            end
                                            counter_inc = counter_inc + 1;
                                            t_vals = [test_dur(test_row) + counter_inc*int, 200];
                                            target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                            if abs(target_diff(1)) < O2_threshold
                                                target_reached = 1; % End while loop
                                            else
                                                target_reached = 0; % Continue while loop
                                            end
                                            time_actual = test_dur(test_row) + counter_inc*int;
                                            if time_actual > ul_time + 10
                                                break % Exit while loop if run duration exceeds time threshold by 10 s
                                            else
                                            end
                                        end
                                    else % Actual value is higher than target
                                        target_reached = 0; % Initialize to enter while loop
                                        counter_dec = 0;
                                        while target_reached == 0
                                            if target_diff(1) < 5 || time_actual < 5 % Smaller step size if difference is within defined value or adjusted time is less than 5 seconds
                                                int = 1;
                                            else
                                                int = 5;
                                            end
                                            counter_dec = counter_dec + 1;
                                            t_vals = [test_dur(test_row) - counter_dec*int, 200];
                                            target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                            if abs(target_diff(1)) < O2_threshold
                                                target_reached = 1; % End while loop
                                            else
                                                target_reached = 0; % Continue while loop
                                            end
                                            time_actual = test_dur(test_row) - counter_dec*int;
                                            if time_actual < ll_time
                                                break;
                                            else
                                            end
                                        end
                                    end
                                elseif abs(target_diff(1)) < O2_threshold && target_diff(2) == 1 % Difference is within oxygen concentration threshold and inflection point is detected
                                    time_actual = test_dur(test_row);
                                else
                                end
                            else % Solution already found
                            end
                        else % Decreasing oxygen segment
                            while double(isempty(sol_row)) == 1 % Run through possible solution sets until defined criteria is met
                                test_row = find(time_val == time_psol_sort(counter_loop));
                                test_row = test_row(1);
                                test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_inc];
                                pred_1 = predict(model_dec,test_x_1); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_inc];
                                pred_2 = predict(model_dec,test_x_2); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment

                                if target_val(test_row) > pred_2 - O2_threshold && target_val(test_row) < pred_1 + O2_threshold
                                    test_x_3 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), mid_inc];
                                    pred_3 = predict(model_dec,test_x_3); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                    pred_range = [pred_2, pred_3, pred_1]; % GPE prediction range
                                    [pred_sorted, idx] = sort(pred_range);
                                    test_vals = [ul_inc, mid_inc, ll_inc];
                                    f_sorted = test_vals(idx);
                                    if length(unique(round(pred_range,3))) < length(pred_range) || range(pred_range) < range_threshold % There does not seem to be significant influence of inlet oxygen fraction for next segment
                                        sol_f_next = ul_inc;
                                    else
                                        sol_f_next = interp1(pred_sorted,f_sorted,target_val(test_row),'linear','extrap'); % Determine inlet oxygen fraction for next segment using GPE results
                                        if sol_f_next > ul_inc
                                            sol_f_next = ul_inc;
                                        elseif sol_f_next < ll_inc
                                            sol_f_next = ll_inc;
                                        else
                                        end
                                    end

                                    % Adjust time for current segment if necessary to keep segment end target and actual value difference within oxygen concentration threshold
                                    t_vals = [test_dur(test_row), 200]; % Run past current segment to fully capture remnant decrease
                                    f_vals = [sol_f(i), sol_f_next];
                                    toggle = 1; % Define function output to be for decreasing oxygen during segment
                                    target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                    if abs(target_diff(1)) > O2_threshold || target_diff(2) ~= 1  % Difference between segment end target and actual value is not within oxygen concencentration threshold or there is no inflection point
                                        if target_diff(1) > 0 || target_diff(2) ~= 1 % Actual value is higher than target or there is no inflection --> increase run duration
                                            target_reached = 0; % Initialize to enter while loop
                                            counter_inc = 0;
                                            while target_reached == 0 % Smaller step size if difference is within defined value
                                                if target_diff(1) < 5
                                                    int = 1;
                                                else
                                                    int = 5;
                                                end
                                                counter_inc = counter_inc + 1;
                                                t_vals = [test_dur(test_row) + counter_inc*int, 200];
                                                target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                                if abs(target_diff(1)) < O2_threshold
                                                    target_reached = 1; % End while loop
                                                else
                                                    target_reached = 0; % Continue while loop
                                                end
                                                time_actual = test_dur(test_row) + counter_inc*int;
                                                if time_actual > ul_time + 10
                                                    break % Exit while loop if run duration exceeds time threshold by 10 s
                                                else
                                                end
                                            end
                                        else % Actual value is lower than target
                                            target_reached = 0; % Initialize to enter while loop
                                            counter_dec = 0;
                                            while target_reached == 0
                                                if target_diff(1) < 5 || time_actual < 5 % Smaller step size if difference is within defined value or adjusted time is less than 5 seconds
                                                    int = 1;
                                                else
                                                    int = 5;
                                                end
                                                counter_dec = counter_dec + 1;
                                                t_vals = [test_dur(test_row) - counter_dec*int, 200];
                                                target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                                if abs(target_diff(1)) < O2_threshold
                                                    target_reached = 1; % End while loop
                                                else
                                                    target_reached = 0; % Continue while loop
                                                end
                                                time_actual = test_dur(test_row) - counter_dec*int;
                                                if time_actual < ll_time
                                                    break;
                                                else
                                                end
                                            end
                                        end
                                    elseif abs(target_diff(1)) < O2_threshold && target_diff(2) == 1 % Difference is within oxygen concentration threshold and inflection point is detected
                                        time_actual = test_dur(test_row);
                                    else
                                    end
                                    sol_row = 1; % Solution identified; break out of while loop
                                else % No solution found from GPE range comparison; continue with while loop
                                    sol_row = [];
                                    counter_loop = counter_loop + 1;
                                    if counter_loop > length(time_psol_sort)
                                        break % No solution found; exit while loop
                                    else
                                    end
                                end
                            end

                            % No solution found from GPE comparison
                            if counter_loop > length(time_psol_sort) % Cycle through all possible solutions in filtered set
                                clear lms
                                for n=1:length(time_psol_sort)
                                    test_row = find(time_val == time_psol_sort(n));
                                    test_row = test_row(1);
                                    test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_inc];
                                    pred_1 = predict(model_inc,test_x_1); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                    test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_inc];
                                    pred_2 = predict(model_inc,test_x_2); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                    lms(n) = (target_val(test_row) - pred_1)^2 + (target_val(test_row) - pred_2)^2;
                                end
                                min_lms = find(lms == min(lms)); % Find solution set with lowest difference between required and predicted remnant increase
                                test_row = find(time_val == time_psol_sort(min_lms));
                                test_row = test_row(1);
                                test_x_1 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ll_inc];
                                pred_1 = predict(model_inc,test_x_1); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                test_x_2 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), ul_inc];
                                pred_2 = predict(model_inc,test_x_2); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                test_x_3 = [ic.c3_20_ic*1e9, sol_f(i), test_dur(test_row), mid_inc];
                                pred_3 = predict(model_inc,test_x_3); % Predict remnant decrease using trained GPE, with assumption for inlet oxygen fraction of the next segment
                                pred_range = [pred_2, pred_3, pred_1]; % GPE prediction range
                                [pred_sorted, idx] = sort(pred_range);
                                test_vals = [ul_inc, mid_inc, ll_inc];
                                f_sorted = test_vals(idx);
                                if length(unique(round(pred_range,3))) < length(pred_range) || range(pred_range) < range_threshold % There does not seem to be significant influence of inlet oxygen fraction for next segment
                                    sol_f_next = ul_inc;
                                else
                                    sol_f_next = interp1(pred_sorted,f_sorted,target_val(test_row),'linear','extrap'); % Determine inlet oxygen fraction for next segment using GPE results
                                    if sol_f_next > ul_inc
                                        sol_f_next = ul_inc;
                                    elseif sol_f_next < ll_inc
                                        sol_f_next = ll_inc;
                                    else
                                    end
                                end

                                % Adjust time for current segment if necessary to keep segment end target and actual value difference within oxygen concentration threshold
                                t_vals = [test_dur(test_row), 200]; % Run past current segment to fully capture remnant increase
                                f_vals = [sol_f(i), sol_f_next];
                                toggle = 1; % Define function output to be for decreasing segment
                                target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                if abs(target_diff(1)) > O2_threshold || target_diff(2) ~= 1 % Difference between segment end target and actual value is not within oxygen concentration threshold or there is no inflection point
                                    if target_diff(1) > 0 || target_diff(2) ~= 1
                                        target_reached = 0; % Initialize to enter while loop
                                        counter_inc = 0;
                                        while target_reached == 0
                                            if target_diff(1) < 5 % Smaller step size if difference is within defined value
                                                int = 1;
                                            else
                                                int = 5;
                                            end
                                            counter_inc = counter_inc + 1;
                                            t_vals = [test_dur(test_row) + counter_inc*int, 200];
                                            target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                            if abs(target_diff(1)) < O2_threshold
                                                target_reached = 1; % End while loop
                                            else
                                                target_reached = 0; % Continue while loop
                                            end
                                            time_actual = test_dur(test_row) + counter_inc*int;
                                            if time_actual > ul_time + 10
                                                break % Exit while loop if run duration exceeds time threshold by 10 s
                                            else
                                            end
                                        end
                                    else % Actual value is lower than target
                                        target_reached = 0; % Initialize to enter while loop
                                        counter_dec = 0;
                                        while target_reached == 0
                                            if target_diff(1) < 5 || time_actual < 5 % Smaller step size if difference is within defined value or adjusted time is less than 5 seconds
                                                int = 1;
                                            else
                                                int = 5;
                                            end
                                            counter_dec = counter_dec + 1;
                                            t_vals = [test_dur(test_row) - counter_dec*int, 200];
                                            target_diff = target_reached_function(k1,k2,k3,k4,k5,k6,Inlet_l2,Inlet_l3,toggle,ic,target,t_vals,f_vals,T_in,target_t);
                                            if abs(target_diff(1)) < O2_threshold
                                                target_reached = 1; % End while loop
                                            else
                                                target_reached = 0; % Continue while loop
                                            end
                                            time_actual = test_dur(test_row) - counter_dec*int;
                                            if time_actual < ll_time
                                                break;
                                            else
                                            end
                                        end
                                    end
                                elseif abs(target_diff(1)) < O2_threshold && target_diff(2) == 1 % Difference is within oxygen concentration threshold and inflection point is detected
                                    time_actual = test_dur(test_row);
                                else
                                end
                            else % Solution already found
                            end
                        end
                    end
                    sol_dur(i) = time_actual;
                end

            else % Run duration is within defined range
            end
        else % Solution determined from primary method
            disp('Primary method successful')
            sol_dur(i) = test_dur(sol_row);
            sol_f_next = test_next(sol_row);
        end
    end
    disp(['iteration = ' num2str(i)])
end

