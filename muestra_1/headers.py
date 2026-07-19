import numpy as np
from astropy.io import fits
import os

def extraer_metadatos_galaxia(nombre_archivo):
    """
    Busca el Centro (xc, yc) y el Position Angle (PA) en el header de un FITS.
    """
    if not os.path.exists(nombre_archivo):
        print(f"Error: El archivo {nombre_archivo} no existe.")
        return None

    with fits.open(nombre_archivo) as hdul:
        # Usualmente los mapas cinemáticos están en la extensión 1 o nombres como 'V_0'
        # Intentamos obtener el header del mapa de velocidad o el primario
        if "V_0" in hdul:
            header = hdul["V_0"].header
            data_shape = hdul["V_0"].data.shape
        else:
            header = hdul[0].header
            data_shape = hdul[0].data.shape

        print(f"\n--- Analizando: {nombre_archivo} ---")

        # 1. BÚSQUEDA DEL CENTRO (En píxeles)
        # CRPIX son las palabras clave estándar de la WCS para el centro
        xc = header.get('CRPIX1')
        yc = header.get('CRPIX2')

        # Si no existen, usamos el centro geométrico de la matriz de datos
        if xc is None or yc is None:
            yc_geo, xc_geo = np.array(data_shape) / 2
            print(f"[*] CRPIX no encontrados. Usando centro geométrico: ({xc_geo}, {yc_geo})")
            xc, yc = xc_geo, yc_geo
        else:
            print(f"[+] Centro encontrado (CRPIX): ({xc}, {yc})")

        # 2. BÚSQUEDA DEL PA (Position Angle)
        # Intentamos varios nombres comunes en archivos de CALIFA/PYCASSO
        posibles_pa = ['PA', 'POSANG', 'PA_OBS', 'PA_OUT', 'PA_PHO', 'PA_FIT']
        pa_encontrado = None

        for clave in posibles_pa:
            if clave in header:
                pa_encontrado = header[clave]
                print(f"[+] PA encontrado en header ('{clave}'): {pa_encontrado}°")
                break
        
        if pa_encontrado is None:
            print("[!] ADVERTENCIA: No se encontró PA en el header.")
            print("    Deberás obtenerlo de la tabla general (califa_master.csv).")

        return {'xc': xc, 'yc': yc, 'pa': pa_encontrado}

# --- EJEMPLO DE USO ---
# Sustituye con el nombre de uno de tus archivos descargados
archivo_test = "K0001_gsd6e.fits"

metadatos = extraer_metadatos_galaxia(archivo_test)

if metadatos:
    print("\nResumen para deproyección:")
    print(f"Centro X: {metadatos['xc']}")
    print(f"Centro Y: {metadatos['yc']}")
    print(f"PA: {metadatos['pa']}")
