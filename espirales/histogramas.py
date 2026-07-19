import pandas as pd
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt

# CARGA Y PREPARACIÓN DE DATOS
df_gme = Table.read("gme.fits").to_pandas()
df_cbe = Table.read("cbe.fits").to_pandas()
df_geo = pd.read_csv("growth_curves.csv", comment='#', header=None).rename(columns={0: 'califaID', 2: 'ba'})
df_morf = pd.read_csv("morphology.csv", comment='#', header=None).rename(columns={0: 'califaID', 5: 'hubtyp', 6: 'subtyp'})

df_gme['califaID'] = df_gme['califaID'].astype(int)
df_cbe['califaID'] = df_cbe['califaID'].astype(int)

# Unión maestra de catálogos
df_pyc = pd.merge(df_gme, df_cbe, on='califaID', suffixes=('_gme', '_cbe'))
df_total = pd.merge(df_pyc, df_morf[['califaID', 'hubtyp', 'subtyp']], on='califaID')
df_total = pd.merge(df_total, df_geo[['califaID', 'ba']], on='califaID')

# Parámetros físicos
q0 = 0.0
mask_sigma = (df_total['sigma_V_gme'] > 80) & (df_total['sigma_V_cbe'] > 80)

def calcular_i(ba_series):
    cos_i2 = (ba_series**2 - q0**2) / (1 - q0**2)
    return np.degrees(np.arccos(np.sqrt(cos_i2.clip(lower=0))))

# Muestra 1: Todas las marcadas como 'S' (incluye potenciales S0)
mask_con_s0 = df_total['hubtyp'].astype(str).str.contains('S', na=False)
df_all_s = df_total[mask_con_s0].copy()
df_filt_s = df_total[mask_con_s0 & mask_sigma].copy()

# Muestra 2: Espirales puras (sin '0' en el subtipo)
mask_sin_s0 = mask_con_s0 & (~df_total['subtyp'].astype(str).str.contains('0', na=False))
df_all_pure = df_total[mask_sin_s0].copy()
df_filt_pure = df_total[mask_sin_s0 & mask_sigma].copy()

# 3. CREACIÓN DEL GRÁFICO MAESTRO
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

def plot_hist(ax, data_all, data_filt, title):
    i_all = calcular_i(data_all['ba'])
    i_filt = calcular_i(data_filt['ba'])
    
    ax.hist(i_all, bins=20, color='lightgray', alpha=0.5, label=f'Total (N={len(i_all)})')
    ax.hist(i_filt, bins=20, color='#2c3e50', alpha=0.8, label=rf'$\sigma > 80$ km/s (N={len(i_filt)})')
    ax.axvline(i_all.mean(), color='gray', linestyle='--', label=f'Media: {i_all.mean():.1f}°')
    ax.axvline(i_filt.mean(), color='#e74c3c', linestyle='-', label=f'Media Filt: {i_filt.mean():.1f}°')
    
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.set_xlabel("Inclinación $i$ [grados]")
    ax.legend(fontsize=9)
    ax.grid(axis='y', alpha=0.2)

# Panel Izquierdo: Con S0
plot_hist(ax1, df_all_s, df_filt_s, "Muestra A: Incluyendo Galaxias S0")
ax1.set_ylabel("Número de Galaxias")

# Panel Derecho: Sin S0
plot_hist(ax2, df_all_pure, df_filt_pure, "Muestra B: Solo Espirales Puras (Sa-Sd)")

plt.suptitle("Comparativa de Selección de Muestra: Impacto de S0 y Límite Cinemático", fontsize=15, y=1.05)
plt.tight_layout()
plt.savefig("histogramas_S.png", dpi=300, bbox_inches='tight')
plt.show()

print(f"Diferencia de galaxias totales: {len(df_all_s) - len(df_all_pure)} lenticulares eliminadas.")
