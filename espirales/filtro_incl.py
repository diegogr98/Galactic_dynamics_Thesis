import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.table import Table
import os

# 1. CARGA DE CATÁLOGOS
df_gme = Table.read("gme.fits").to_pandas()
df_cbe = Table.read("cbe.fits").to_pandas()
df_geo = pd.read_csv("growth_curves.csv", comment='#', header=None).rename(columns={0: 'califaID', 1: 'pa', 2: 'ba'})
df_morf = pd.read_csv("morphology.csv", comment='#', header=None).rename(columns={0: 'califaID', 5: 'hubtyp', 6: 'subtyp'})

# Limpieza y unión
df_gme['califaID'] = df_gme['califaID'].astype(int)
df_cbe['califaID'] = df_cbe['califaID'].astype(int)

df_pyc = pd.merge(df_gme, df_cbe, on='califaID', suffixes=('_gme', '_cbe'))
df_total = pd.merge(df_pyc, df_morf[['califaID', 'hubtyp', 'subtyp']], on='califaID')
df_total = pd.merge(df_total, df_geo[['califaID', 'pa', 'ba']], on='califaID')

# calculando inclinaciones
q0 = 0.0 
cos_i2 = (df_total['ba']**2 - q0**2) / (1 - q0**2)
df_total['inclination'] = np.degrees(np.arccos(np.sqrt(cos_i2.clip(lower=0))))

# Máscara estricta: i < 30, Sigma > 80
mask_sigma = (df_total['sigma_V_gme'] > 80)
mask_inc = (df_total['inclination'] < 30)
mask_espiral = df_total['hubtyp'].astype(str).str.contains('S', na=False)
mask_sin_s0 = ~df_total['subtyp'].astype(str).str.contains('0', na=False)

df_espirales = df_total[mask_sigma & mask_inc & mask_espiral & mask_sin_s0].copy()

# 4. GUARDAR CSV FINAL
columnas_finales = ['califaID', 'hubtyp', 'subtyp', 'inclination', 'pa', 'sigma_V_gme']
df_espirales[columnas_finales].to_csv("espirales_final.csv", index=False)

print(f"\nÉxito. Se encontraron {len(df_espirales)} espirales puras.")
print("Archivo guardado como: espirales_final.csv")
