import csv
import math
import os
import numpy as np

repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
csv_path = os.path.join(repo_root, 'Compiled data', 'revision_tutor', 'GSD_sensitivity_summary.csv')

existing = []
with open(csv_path, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        existing.append(row)

gsd_known = np.array([60, 80, 100])
seca_known = np.array([3.3417, 1.3369, 0.7410])
comb_known = np.array([0.8714, 0.4914, 0.3447])

log_seca = np.log(seca_known)
log_comb = np.log(comb_known)
B, log_A = np.polyfit(log_seca, log_comb, 1)
A = math.exp(log_A)
print(f"Fuel fit: comb = {A:.6f} * seca^{B:.6f}")

rows = []
for gsd in range(30, 121, 10):
    D_mm = round(2273 / gsd)
    masa_seca = 3.3417 * (60 / gsd) ** 3
    masa_comb = A * (masa_seca ** B)
    masa_total_sat = masa_seca + masa_comb
    masa_total_const = 2 * masa_total_sat  # 2 sat
    detectores = math.ceil(170000 / (gsd * 2048))
    
    rows.append({
        'GSD_m': str(gsd),
        'Num_Satelites': '2',
        'Num_Telescopios': '2',
        'Detectores_Necesarios': str(detectores),
        'Altura_km': '520',
        'Diametro_Pupila_mm': str(D_mm),
        'Swath_Requerido_km': '170',
        'Masa_Seca_Satelite_kg': f'{masa_seca:.6f}',
        'Masa_Combustible_kg': f'{masa_comb:.6f}',
        'Masa_Total_Constelacion_kg': f'{masa_total_const:.6f}',
        'Tipo_Telescopio': 'Refractivo',
    })

fieldnames = ['GSD_m', 'Num_Satelites', 'Num_Telescopios', 'Detectores_Necesarios',
              'Altura_km', 'Diametro_Pupila_mm', 'Swath_Requerido_km',
              'Masa_Seca_Satelite_kg', 'Masa_Combustible_kg', 'Masa_Total_Constelacion_kg',
              'Tipo_Telescopio']

with open(csv_path, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} rows to {csv_path}")
for r in rows:
    print(f"  GSD={r['GSD_m']:>3}m  D={r['Diametro_Pupila_mm']:>3}mm  seca={float(r['Masa_Seca_Satelite_kg']):.4f}kg  comb={float(r['Masa_Combustible_kg']):.4f}kg  total_const={float(r['Masa_Total_Constelacion_kg']):.4f}kg")
