import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

data = pd.read_csv('results/heatmap_data.csv')

x_coords = data['μϵ_W'].values
y_coords = data['ΔV'].values
values = data['P_w2'].values
#values = data['η_ηc'].values

# Create a meshgrid for the coordinates
X, Y = np.meshgrid(np.unique(x_coords), np.unique(y_coords))
Z = values.reshape(X.shape)

# Custom colormap
cmap = LinearSegmentedColormap.from_list(
    'custom_cmap', [(0, 'white'), (0.1, 'white'), (0.25, 'blue'), (0.5, 'yellow'), (0.75, 'red'), (1.0, 'black')])

cmap.set_bad('white')

contour = plt.contourf(X, Y, Z, cmap=cmap)
plt.colorbar(label=r'$P/W^2$')
plt.xlabel(r'Av. chemical potential $(\mu - \epsilon)/W$')
plt.ylabel(r'Potential difference $\Delta V/W$')
plt.show()
