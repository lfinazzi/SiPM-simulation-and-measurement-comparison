#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 25 13:06:20 2025

@author: finazzi
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson

def sipm_pulse(t, amplitude=1.0, rise_time=5, fall_time=30, pulse_width=100, baseline=0):
    """
    Simulate an SiPM pulse with exponential rise and fall.

    Parameters:
        t (array): Time array in ns
        amplitude (float): Peak amplitude of the pulse
        rise_time (float): Time for exponential rise in ns
        fall_time (float): Time for exponential fall in ns
        pulse_width (float): Total width (FWHM-like) of the pulse in ns
        baseline (float): Baseline offset

    Returns:
        pulse (array): SiPM pulse signal
    """
    center = pulse_width / 2
    pulse = amplitude * (1 - np.exp(-(t - center) / rise_time)) * np.exp(-(t - center) / fall_time)
    pulse[t < center] = 0  # Zero out non-causal part
    return pulse + baseline

def plot_and_integrate_pulse(amplitude=1.0, rise_time=5, fall_time=30,
                             pulse_width=100, baseline=0, time_range=(0, 200), dt=0.5):
    """
    Generate, plot, and integrate an SiPM pulse.

    Parameters:
        amplitude, rise_time, fall_time, pulse_width, baseline: SiPM parameters
        time_range (tuple): Start and end of time window in ns
        dt (float): Time resolution in ns
    """
    t = np.arange(time_range[0], time_range[1], dt)
    pulse = sipm_pulse(t, amplitude, rise_time, fall_time, pulse_width, baseline)

    # Integrate pulse (area under the curve - baseline corrected)
    charge = simpson(pulse - baseline, t)  # ns·units (e.g., ns·mV)

    # Plotting
    plt.figure(figsize=(10, 4))
    plt.plot(t, pulse, label='SiPM Pulse', color='blue')
    plt.axhline(y=baseline, color='gray', linestyle='--', label='Baseline')
    plt.title('Simulated SiPM Pulse')
    plt.xlabel('Time (ns)')
    plt.ylabel('Amplitude (a.u.)')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    print(f"Integrated charge: {charge:.3f} ns·a.u.")

# Example usage
plot_and_integrate_pulse(
    amplitude=1.15,
    rise_time=15,
    fall_time=30,
    pulse_width=150,
    baseline=0.05,
    time_range=(0, 200),
    dt=0.5
)
