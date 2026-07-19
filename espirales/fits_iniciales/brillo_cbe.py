import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import os

# --- CONFIGURACIÓN GENERAL ---
RUTA_CSV_COMPLETO = os.path.expanduser("~/califa/espirales/espirales_final_con_bulbo.csv")
TAMANO_BIN = 2.0    # Tamaño del bin en píxeles
FILTRO_SN = 5       # Umbral de Señal/Ruido para limpiar los bordes externos

# Nombre de la extensión exacta para este análisis
EXT_FITS = "L_5635"
# Corregido: Removido \text{} para evitar el fallo del renderizador de Matplotlib
Y_LABEL = r"$\log_{10}(L_{5635}\AA)$"

# 1. CARGAR DATOS DESDE EL CSV UNIFICADO
df_galaxias = pd.read_csv(RUTA_CSV_COMPLETO)
df_galaxias.columns = df_galaxias.columns.str.strip()

lista_ids = df_galaxias['califaID'].tolist()
print(f"Procesando {len(lista_ids)} galaxias para el perfil de: {EXT_FITS}")

# --- MÉTODOS DE EXTRACCIÓN DE RADIOS ---

def calcular_radios(archivo, pa_deg, inc_deg, ext_nombre):
    if not os.path.exists(archivo):
        return None, None, None

    phi = np.radians(pa_deg + 90)
    cos_i = np.cos(np.radians(inc_deg))

    with fits.open(archivo) as hdul:
        if ext_nombre not in hdul:
            print(f"Advertencia: No se encontró la extensión {ext_nombre} en {archivo}")
            return None, None, None
            
        data_map = hdul[ext_nombre].data.copy()
        sn = hdul["SN"].data if "SN" in hdul else None
        header = hdul[ext_nombre].header

        # Máscara usando el canal de Señal/Ruido de tus cubos
        if sn is not None:
            data_map[sn < FILTRO_SN] = np.nan

        # Centro astronómico (CRPIX) tomado del header del FITS
        xc = header.get('CRPIX1', data_map.shape[1] / 2)
        yc = header.get('CRPIX2', data_map.shape[0] / 2)

        # Crear mallas de coordenadas espaciales
        y, x = np.indices(data_map.shape)
        xt, yt = x - xc, y - yc

        # Modelado de geometrías radiales (Circular vs Elíptica Desinclinada)
        r_circular = np.sqrt(xt**2 + yt**2)
        x_prime = xt * np.cos(phi) + yt * np.sin(phi)
        y_prime = -xt * np.sin(phi) + yt * np.cos(phi)
        r_elip_desinc = np.sqrt(x_prime**2 + (y_prime / cos_i)**2)

        return r_circular.flatten(), r_elip_desinc.flatten(), data_map.flatten()

def binning_radial(r, data, step):
    bins = np.arange(0, np.nanmax(r), step)
    bc, bv = [], []
    for i in range(len(bins)-1):
        mask = (r >= bins[i]) & (r < bins[i+1])
        if np.any(mask) and not np.all(np.isnan(data[mask])):
            # Mediana robusta en cada anillo de píxeles
            mediana = np.nanmedian(data[mask])
            if not np.isnan(mediana):
                bc.append((bins[i] + bins[i+1]) / 2)
                bv.append(mediana)
    return np.array(bc)[:len(bv)], np.array(bv)


# --- CREACIÓN DEL MOSAICO DE GRÁFICAS ---
num_galaxias = len(lista_ids)
columnas = 3
filas = int(np.ceil(num_galaxias / columnas))

fig, axes = plt.subplots(filas, columnas, figsize=(12, 3 * filas), sharex=True)
axes = axes.flatten()

for i, id_obj in enumerate(lista_ids):
    id_str = f"{id_obj:04d}"
    archivo_fits = f"K{id_str}_zca6e.fits"

    fila_galaxia = df_galaxias[df_galaxias['califaID'] == id_obj]
    pa = fila_galaxia['pa'].values[0]
    inc = fila_galaxia['inclination'].values[0]
    re_g = fila_galaxia['radio_bulbo_g'].values[0]
    re_i = fila_galaxia['radio_bulbo_i'].values[0]

    ax = axes[i]

    # Extraer perfiles utilizando LOG_Z_FLUX
    r_circ, r_el_des, y_data = calcular_radios(archivo_fits, pa, inc, EXT_FITS)

    if r_circ is not None:
        bc_circ, y_circ = binning_radial(r_circ, y_data, step=TAMANO_BIN)
        bc_el_des, y_el_des = binning_radial(r_el_des, y_data, step=TAMANO_BIN)

        # Graficar curvas radiales (ya están en logaritmo directo)
        if len(y_circ) > 0:
            ax.plot(bc_circ, y_circ, 'o-', color='darkcyan', markersize=3, linewidth=1, label='cbe Circular')
        if len(y_el_des) > 0:
            ax.plot(bc_el_des, y_el_des, 's-', color='indigo', markersize=3, linewidth=1, label='cbe Elíptico')

        # --- LÍNEAS VERTICALES DE LOS RADIOS DEL BULBO (BANDA g e i) ---
        if re_g is not None and not np.isnan(re_g):
            ax.axvline(re_g, color='crimson', linestyle=':', linewidth=1.5, label=f'$R_e$ g ({re_g:.1f}″)')
        if re_i is not None and not np.isnan(re_i):
            ax.axvline(re_i, color='darkorange', linestyle='--', linewidth=1.2, label=f'$R_e$ i ({re_i:.1f}″)')

        ax.set_xlim(0, 40)
        ax.set_title(f"K{id_str} ($i$={int(inc)}°)", fontsize=10)
        ax.grid(True, alpha=0.2)

        if i >= (filas - 1) * columnas: ax.set_xlabel('Radio [píxeles]')
        if i % columnas == 0: ax.set_ylabel(Y_LABEL)
        ax.legend(loc='upper right', fontsize=5.5, framealpha=0.6)

    else:
        ax.text(0.5, 0.5, f"FITS K{id_str}\nsin {EXT_FITS}", ha='center', va='center')
        ax.set_title(f"K{id_str} (Faltante)")

for j in range(num_galaxias, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
nombre_salida = 'cbe_flux.png'
plt.savefig(nombre_salida, dpi=200)
print(f"¡Proceso terminado con éxito! Mosaico guardado como '{nombre_salida}'")
plt.show()
