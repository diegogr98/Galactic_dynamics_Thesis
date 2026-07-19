echo "Descargando archivos necesarios para el análisis..."

# Archivos de Pycasso
wget -O gme.fits https://pycasso.iaa.es/tables/pycasso_integrated_main_gsd6e.fits
wget -O cbe.fits https://pycasso.iaa.es/tables/pycasso_integrated_main_zca6e.fits

# Archivos de Califa
wget -O growth_curves.csv https://califa.caha.es/FTP-PUB/docs/DR3/CALIFA_4_MS_GC_mag.csv
wget -O morphology.csv https://califa.caha.es/FTP-PUB/docs/DR3/CALIFA_2_MS_class.csv

echo "---Archivos descargados correctamente---"
