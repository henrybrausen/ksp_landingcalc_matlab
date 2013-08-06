ksp_landingcalc_matlab
======================

MATLAB Scripts for Calculating Landing Trajectories in Kerbal Space Program

Calculates what periapsis height will let you land directly under your periapsis for atmospheric bodies.

Edit the value of "dt" in landingcalc3.m to change the numerical integration accuracy.
Run enhancedlandingplots.m to generate the plots and data for each body.

These calculations can take a long time to run, depending on the accuracy and number of data points you want.
On my machine, each data point took on the order of 3 seconds to compute, so keep that in mind!
