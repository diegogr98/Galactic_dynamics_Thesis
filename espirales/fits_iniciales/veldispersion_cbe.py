import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import os

# --- CONFIGURACIÓN GENERAL ---
# Ahora solo necesitamos un único archivo CSV que ya tiene todo integrado
RUTA_CSV_COMPLETO = os.path.expanduser("~/califa/espirales/espirales_final_con_bulbo.csv")
TAMANO_BIN = 2.0    # Tamaño del bin solicitado para todos por igual (2)
FILTRO_SN = 5       # Umbral de Señal/Ruido para enmascarar píxeles malos

# 1. CARGAR TODOS LOS DATOS DESDE TU NUEVO CSV UNIFICADO
df_galaxias = pd.read_csv(RUTA_CSV_COMPLETO)
df_galaxias.columns = df_galaxias.columns.str.strip()

lista_ids = df_galaxias['califaID'].tolist()
print(f"Procesando {len(lista_ids)} galaxias usando el archivo unificado con bandas g e i...")

# --- MÉTODOS DE EXTRACCIÓN DE RADIOS ---

def calcular_radios(archivo, pa_deg, inc_deg):
    if not os.path.exists(archivo):
        return None, None, None

    phi = np.radians(pa_deg + 90)
    cos_i = np.cos(np.radians(inc_deg))

    with fits.open(archivo) as hdul:
        ext = "V_D" if "V_D" in hdul else "SIG_V"
        sigma = hdul[ext].data.copy()
        sn = hdul["SN"].data
        header = hdul[ext].header

        # Máscara de calidad
        sigma[sn < FILTRO_SN] = np.nan

        # Centro astronómico desde las palabras clave del Header
        xc = header.get('CRPIX1', sigma.shape[1] / 2)
        yc = header.get('CRPIX2', sigma.shape[0] / 2)

        # Crear la malla de píxeles
        y, x = np.indices(sigma.shape)
        xt, yt = x - xc, y - yc

        # --- APROXIMACIÓN FACE-ON (Distancia Circular Directa) ---
        r_circular = np.sqrt(xt**2 + yt**2)

        # Rotación para alinear con el Ángulo de Posición (PA)
        x_prime = xt * np.cos(phi) + yt * np.sin(phi)
        y_prime = -xt * np.sin(phi) + yt * np.cos(phi)

        # --- APROXIMACIÓN ELÍPTICA DEPRIMIDA (Desinclinando la galaxia) ---
        r_elip_desinc = np.sqrt(x_prime**2 + (y_prime / cos_i)**2)

        return r_circular.flatten(), r_elip_desinc.flatten(), sigma.flatten()

def binning_radial(r, sig, step):
    bins = np.arange(0, np.nanmax(r), step)
    bc, bv = [], []
    for i in range(len(bins)-1):
        mask = (r >= bins[i]) & (r < bins[i+1])
        if np.any(mask) and not np.all(np.isnan(sig[mask])):
            mediana = np.nanmedian(sig[mask]) # Mediana robusta
            if mediana > 0:
                bc.append((bins[i] + bins[i+1]) / 2)
                bv.append(mediana)
    return np.array(bc)[:len(bv)], np.array(bv)


# --- CREACIÓN DEL MOSAICO DE GRÁFICAS ---
num_galaxias = len(lista_ids)
columnas = 3
filas = int(np.ceil(num_galaxias / columnas))

# Dimensiones del mosaico
fig, axes = plt.subplots(filas, columnas, figsize=(12, 3 * filas), sharex=True)
axes = axes.flatten() # Aplanar para el loop

for i, id_obj in enumerate(lista_ids):
    id_str = f"{id_obj:04d}"
    archivo_fits = f"K{id_str}_zca6e.fits"

    # Extraer TODOS los datos de la fila de esta galaxia directamente
    fila_galaxia = df_galaxias[df_galaxias['califaID'] == id_obj]
    
    pa = fila_galaxia['pa'].values[0]
    inc = fila_galaxia['inclination'].values[0]
    
    # Extraemos los dos radios que ya guardamos previamente en el CSV unificado
    re_g = fila_galaxia['radio_bulbo_g'].values[0]
    re_i = fila_galaxia['radio_bulbo_i'].values[0]

    ax = axes[i]

    # Procesar archivo FITS usando el Header
    r_circ, r_el_des, s_data = calcular_radios(archivo_fits, pa, inc)

    if r_circ is not None:
        bc_circ, bv_circ = binning_radial(r_circ, s_data, step=TAMANO_BIN)
        bc_el_des, bv_el_des = binning_radial(r_el_des, s_data, step=TAMANO_BIN)

        # Graficar perfiles radiales comparativos simultáneos
        if len(bv_circ) > 0:
            ax.plot(bc_circ, np.log10(bv_circ), 'o-', color='darkcyan',
                    markersize=3, linewidth=1, label='cbe Circular')
        if len(bv_el_des) > 0:
            ax.plot(bc_el_des, np.log10(bv_el_des), 's-', color='indigo',
                    markersize=3, linewidth=1, label='cbe Elíptico')

        # Línea horizontal del límite cinemático de referencia
        ax.axhline(np.log10(80), color='teal', linestyle='--', alpha=0.6, label='Límite ~80 km/s')

        # --- LÍNEAS VERTICALES DESDE EL CSV UNIFICADO ---
        # Línea vertical para el Radio en Banda g (Línea punteada carmesí)
        if re_g is not None and not np.isnan(re_g):
            ax.axvline(re_g, color='crimson', linestyle=':', linewidth=1.5,
                       label=f'$R_{{e}}$ g ({re_g:.1f}″)')

        # Línea vertical para el Radio en Banda i (Línea segmentada naranja oscuro)
        if re_i is not None and not np.isnan(re_i):
            ax.axvline(re_i, color='darkorange', linestyle='--', linewidth=1.2,
                       label=f'$R_{{e}}$ i ({re_i:.1f}″)')

        ax.set_ylim(0.5, 2.5)
        ax.set_xlim(0, 40)

        ax.set_title(f"Galaxia K{id_str} ($i$={int(inc)}°)", fontsize=10)
        ax.grid(True, alpha=0.2)

        # Etiquetas de los ejes distribuidas de forma limpia
        if i >= (filas - 1) * columnas: ax.set_xlabel('Radio [píxeles]')
        if i % columnas == 0: ax.set_ylabel(r'Log_10$\sigma$ [km/s]')

        # Leyenda pequeña para que quepa en los subplots
        ax.legend(loc='lower left', fontsize=5.5, framealpha=0.6)

    else:
        ax.text(0.5, 0.5, f"FITS K{id_str}\nno encontrado", ha='center', va='center')
        ax.set_title(f"K{id_str} (Faltante)")

# Apagar subplots vacíos si corresponde
for j in range(num_galaxias, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.savefig('vel_disp_cbe.png', dpi=200)
print("¡Mosaico final generado con éxito sin usar archivos externos de fotometría!")
plt.show()
