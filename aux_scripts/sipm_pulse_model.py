import numpy as np
import matplotlib.pyplot as plt

# Pulse parameters
t0 = 10e-9       # Pulse start time [s]
tau_rise = 20e-9  # Rise time constant [s]
tau_fall = 40e-9 # Fall time constant [s]

# Normalized SiPM pulse (area = 1)
def sipm_pulse(t, t0, tau_rise, tau_fall):
    dt = t - t0
    pulse = np.where(t < t0, 0.0, (np.exp(-dt / tau_fall) - np.exp(-dt / tau_rise)) / (tau_fall - tau_rise))
    return pulse

# Cumulative integral of the pulse from t0 to t
def sipm_pulse_integral(t, t0, tau_rise, tau_fall):
    dt = t - t0
    integral = np.where(
        t < t0,
        0.0,
        (tau_fall * (1 - np.exp(-dt / tau_fall)) - tau_rise * (1 - np.exp(-dt / tau_rise))) / (tau_fall - tau_rise)
    )
    return integral

# Time vector
t = np.linspace(0, 500e-9, 1000)

# Evaluate pulse and integral
v = sipm_pulse(t, t0, tau_rise, tau_fall)
integ = sipm_pulse_integral(t, t0, tau_rise, tau_fall)

# Plotting
plt.figure(figsize=(10, 5))

plt.subplot(1, 2, 1)
plt.plot(t * 1e9, v, label='Normalized pulse')
plt.xlabel('Time [ns]')
plt.ylabel('Amplitude [1/s]')
plt.title('SiPM Pulse Shape')
plt.grid(True)
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(t * 1e9, integ, label='Cumulative integral', color='orange')
plt.xlabel('Time [ns]')
plt.ylabel('Integrated area')
plt.title('Pulse Integral')
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.show()
