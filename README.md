# Predictive-Corrective-Wave-Forecasting

## Content 
Contains the codes to carry out corrective wave forecasting using ocean wave buoy data. Uses *Conservative Corrected Smoothed Particle Hydrodynamics (CCSPH)* to simulate ocean waves based on ocean wave spectra obtained from wave buoys at different locations. Following this, the model simulation results and observations at the wave buoys are used to train a range of neural networks to carry out wave elevation forecasting over a certain horizon ranging from 5 minutes to 1 hour. 
To investigate the effect of data assimilation in this process, data assimilation techniques including optimal interpolation, perturbed observation ensemble kalaman filter (PO-EnKF), local ensemble transform kalman filter (LETKF) have been also used using the observations and the neural network predictions at different intervals durinf the entire forecasting horizon. This provides the corrected predictions which can be used for the predictions over the next interval.
The facilitative impact of this corrective approach could be examined in these studies, with the approach achieving 25% reduction in MAE compared to the sole forecasting process.
