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
  1. Try other minimization scalars -> Currently, NLL was miscalculated. Fix this for improved minimization. The code was using p*log(k) instead of k*log(p) in NLL, 
  2. Implement a second minimization using the timing histogram to fine tune timing parameters. First, minimization of certain parameters will be performed using the finger spectrum (with timing parameters fixed). Then, the timing histogram will be used to fine tune timing parameters, maintaining all other parameters fixed.
  3. Implement a temporal signal/waveform generation pipeline (can be constructed with a pulse template of desired SiPM),
  4. Figure out what causes shoulders in finger spectrum. Simulation seems to indicate that it's caused by afterpulsing, but I'm not sure.

