from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

hdul = fits.open("/home/diego/califa/muestra_espirales/fits_iniciales/K0311_gsd6e.fits")

sigma = hdul["V_D"].data
sn = hdul["SN"].data

# máscara de calidad
sigma[sn < 5] = np.nan

plt.imshow(sigma, origin='lower', cmap='inferno')
plt.colorbar(label='Dispersión (km/s)')
plt.title('Mapa de dispersión de velocidades')

plt.show()
