To run gas phase code with a sample parameter set:
1. Download:
   - chamber_gas_phase.m
   - Chamber_gas_parameters_filter_PAPER.m (from Additional Codes folder)
   - x_pars_post_gas_final.csv (from Additional Codes folder)
   - y_out_post_gas_final.csv (from Additional Codes folder)
   - Flow_data.mat (from Model Calibration folder)
2. Run Chamber_gas_parameters_filter_PAPER.m
3. Run chamber_gas_phase.m 

To run gas-liquid phase code with a sample parameter set:
1. Download:
   - chamber_gas_and_liquid_phase.m
   - Chamber_liquid_parameters_filter_PAPER.m (from Additional Codes folder)
   - x_pars_post_liquid_final.csv (from Additional Codes folder)
   - y_out_post_liquid_final.csv (from Additional Codes folder)
   - Flow_data.mat (from Model Calibration folder)
2. Run Chamber_liquid_parameters_filter_PAPER.m
3. Run chamber_gas_and_liquid_phase.m 

NOTES:
- A sample IH protocol (Sample_IH_Protocol.csv) has been provided for IH run testing
- The codes automatically choose a random filtered parameter set to run the model. This can be adjusted in the first section of both chamber_gas_phase.m and chamber_gas_and_liquid_phase.m
- Parameters such as starting liquid volume, cell density, cellular metabolic coefficients, etc. can be adjusted within the codes
