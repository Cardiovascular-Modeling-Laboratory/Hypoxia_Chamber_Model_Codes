1. Download:
   - Chamber_liquid_parameters_filter_PAPER.m (from Additional Codes folder)
   - x_pars_post_liquid_final.csv (from Additional Codes folder)
   - y_out_post_liquid_final.csv (from Additional Codes folder)
   - Flow_data.mat (from Model Calibration folder)
   - opt_function_primary.m (Optimization Codes folder)
   - opt_function_secondary.m (Optimization Codes folder)
   - target_reached_function.m (Optimization Codes folder)
   - Chamber_initial_condition.m (Optimization Codes folder)
   - Training_Data_Inc_Dec.xlsx (Optimization Codes folder)
2. Run Chamber_liquid_parameters_filter_PAPER.m to obtain a parameter set for optimization
3. Run Protocol_predicition_PAPER.m to run optimization for protocol prediction

NOTES:
- The optimization code automatically chooses a random filtered parameter set to run the model. This can be adjusted in the first section of both Protocol_prediction_PAPER.m
- The file target_pattern_PAPER.m was used and is provided as a sample to run the optimization code 
