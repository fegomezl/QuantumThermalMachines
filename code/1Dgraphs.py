#!/usr/bin/python3
import pandas as pd
import matplotlib.pyplot as plt

Temp50 = pd.read_csv('results/temp_50.csv')
Temp100 = pd.read_csv('results/temp_100.csv')
Temp200 = pd.read_csv('results/temp_200.csv')

Gamma50 = pd.read_csv('results/gamma_50.csv')
Gamma100 = pd.read_csv('results/gamma_100.csv')
Gamma200 = pd.read_csv('results/gamma_200.csv')

Energy20 = pd.read_csv('results/energy_20.csv')
Energy50 = pd.read_csv('results/energy_50.csv')
Energy100 = pd.read_csv('results/energy_100.csv')

fig, (ax1,ax2,ax3) = plt.subplots(1,3,sharey=True)
ax1.set_ylabel(r'Particle Current $J_p$')
ax2.set_ylabel(r'Particle Current $J_p$')
ax3.set_ylabel(r'Particle Current $J_p$')
ax1.grid()
ax2.grid()
ax3.grid()
ax1.plot(Temp50['T'], Temp50['Jₚ'])
ax1.plot(Temp100['T'], Temp100['Jₚ'])
ax1.plot(Temp200['T'], Temp200['Jₚ'])
plt.show()
