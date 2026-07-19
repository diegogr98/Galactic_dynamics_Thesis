import pandas as pd
from astropy.table import Table

# 1. Cargar las tablas de los dos modelos
try:
    df_gme = Table.read("gme.fits").to_pandas()
    df_cbe = Table.read("cbe.fits").to_pandas()
    
    print("--- COLUMNAS ENCONTRADAS EN GMe ---")
    print(df_gme.columns.tolist())
    
    print("\n--- PRIMERAS FILAS PARA VER VALORES ---")
    print(df_gme.head(20))
    
except Exception as e:
    print(f"Error al leer los archivos: {e}")

