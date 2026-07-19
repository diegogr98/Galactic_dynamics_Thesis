from astropy.io import fits

# Abre el archivo y muestra un resumen de su estructura
with fits.open('K0311_gsd6e.fits') as hdul:
    hdul.info()
