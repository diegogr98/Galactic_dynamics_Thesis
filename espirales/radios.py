import numpy as np
import pandas as pd

print("Cargando espirales_final.csv...")
df_espirales = pd.read_csv("espirales_final.csv")
df_espirales.columns = df_espirales.columns.str.strip()

print("Cargando 2d_descomposicion_fot.csv usando índices posicionales...")
# Cargamos el archivo sin usar nombres de columnas (header=None) ya que no vienen en el formato estándar
df_photdec = pd.read_csv(
    "2d_descomposicion_fot.csv", comment="#", header=None, skipinitialspace=True
)

# Creamos un subconjunto usando los índices numéricos exactos de las columnas:
# Columna 1 (Índice 0)   -> CALIFAID
# Columna 15 (Índice 14) -> Re_g
# Columna 127 (Índice 126) -> Re_i
df_radios = df_photdec[[0, 14, 126]].copy()

# Renombramos las columnas para poder trabajar de forma limpia
df_radios.columns = ["CALIFAID", "radio_bulbo_g", "radio_bulbo_i"]

# Convertimos los IDs a enteros para asegurar un cruce perfecto sin errores de formato
df_espirales["califaID"] = df_espirales["califaID"].astype(int)
df_radios["CALIFAID"] = df_radios["CALIFAID"].astype(int)

print("Efectuando la unión de los catálogos...")
df_final = pd.merge(
    df_espirales,
    df_radios,
    left_on="califaID",
    right_on="CALIFAID",
    how="left",
)

# Eliminamos la columna de ID repetida tras el merge
df_final = df_final.drop(columns=["CALIFAID"])

# Limpiamos los valores nulos astronómicos de no-detección (-999.0)
df_final["radio_bulbo_g"] = df_final["radio_bulbo_g"].replace(-999.0, np.nan)
df_final["radio_bulbo_i"] = df_final["radio_bulbo_i"].replace(-999.0, np.nan)

# Guardar el archivo combinado
df_final.to_csv("espirales_final_con_bulbo.csv", index=False)

print("\n¡Fusión completada exitosamente!")
print("Archivo guardado como: 'espirales_final_con_bulbo.csv'")
print(f"Columnas resultantes: {list(df_final.columns)}")
print("\nMuestra de los primeros datos combinados:")
print(df_final[["califaID", "radio_bulbo_g", "radio_bulbo_i"]].dropna().head())
