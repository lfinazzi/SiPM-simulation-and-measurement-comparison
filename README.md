# SiPM-simulation-and-measurement-comparison
The aim of this repository is to simulate SiPM data and to compare it with real measurements acquired with CAEN DT5751. In addition, a Minimizer was implemented to better tune simulation parameters. Fined tune simulation can then be used for SiPM data generation.

Measurements were made with CAEN DT5751 and these measurements contain charge in a 75 ns gate (among other parameters, like timetag of event or pile up flags).

The simulation part of the code generates SiPM primary, crosstalk and afterpulsing events based on input parameters (like gain, PDE, DCR, crosstalk probability, among others). Each event simulated saves the following data:

  1. Event number,
  2. Timetag,
  3. Charge,
  4. label ("P" for primary, "CT1" for CT of primary, "CT2" for crosstalk of crosstalk or "AP" for afterpulsing).

The code can plot the finger spectrum (using charge values) of the data and simulation in the same plot for direct comparison.

# TODO
  1. Try other minimization scalars (currently using Negative Log Likelihood, but it's erratic),
  2. Refine minimization procedure and take care of noisy gradients,
  3. Implement a way to leave some of the VariableParameters fixed for debug purposes and avoid overfitting,
  4. Implement a temporal signal/waveform generation pipeline (can be constructed with a pulse template of desired SiPM).

