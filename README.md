# SiPM-simulation-and-measurement-comparison
The aim of this repository is to simulate SiPM data and to compare it with real measurements acquired with CAEN DT5751.

Measurements were made with CAEN DT5751 and these measurements contain charge in a 75 ns gate (among other parameters, like timetag of event or pile up flags).

The simulation part of the code generates SiPM primary, crosstalk and afterpulsing events based on input parameters (like gain, PDE, DCR, crosstalk probability, among others). Each event simulated saves the following data:

  1. Event number,
  2. Timetag,
  3. Charge,
  4. label ("P" for primary, "CT1" for CT of primary, "CT2" for crosstalk of crosstalk or "AP" for afterpulsing).

The code can plot the finger spectrum (using charge values) of the data and simulation in the same plot for direct comparison.

# TODO
  1. Implement a Minimizer class that takes a RealData object and a Simulation object and performs minimization on simulation parameters to find the optimal parameters that fit the measured data.
  2. Implement plotting of minimization steps to see how the simulation gets closer to the measured data as iteration of parameter values is performed. 
