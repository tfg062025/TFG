import os
import glob
import numpy as np
import matplotlib.pyplot as plt

carpeta_fitxers = "Fitxers"

noms_parametres = ["x0", "phi", "lambda", "mu", "sigma", "delta"]

carpeta_base = "Histogrames evolució de la població"
os.makedirs(carpeta_base, exist_ok=True)

for nom in noms_parametres:
    os.makedirs(os.path.join(carpeta_base, nom), exist_ok=True)

def extreure_numero(nom_fitxer):
    return int(nom_fitxer.split('_')[-1].split('.')[0])

fitxers = sorted(glob.glob(os.path.join(carpeta_fitxers, "generació_*.txt")), key=extreure_numero)

for idx, fitxer in enumerate(fitxers):
    with open(fitxer, 'r') as f:
        dades = [list(map(float, linia.strip().split())) for linia in f if linia.strip()]
    
    dades = np.array(dades).T 
    generacio_idx = f"{idx:03d}" 

    for i, nom_param in enumerate(noms_parametres):
        plt.figure(figsize=(8, 5))
        plt.hist(dades[i], bins=50, color='skyblue', edgecolor='black')
        plt.title(f"Histograma dels individus de P({generacio_idx})")
        plt.xlabel(nom_param)
        plt.ylabel("Freqüència")
        plt.tight_layout()
        plt.savefig(os.path.join(carpeta_base, nom_param, f"generació_{generacio_idx}.png"))
        plt.close()
