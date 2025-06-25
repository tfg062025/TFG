import numpy as np
import os

def llegir_observacions(nom_fitxer,carpeta):
    observacions = []
    directori = os.path.join(os.path.dirname(__file__), '..', 'Resultats',carpeta, nom_fitxer)
    with open(directori, 'r', encoding='utf-8') as f:
        contingut = f.read()

    blocs = contingut.strip().split('\n\n')
    
    for bloc in blocs:
        files = []
        for linia in bloc.strip().split('\n'):
            if linia.strip():
                nombres = list(map(float, linia.strip().split()))
                files.append(nombres)
        observacions.append(files)

    return observacions

def calcular_mitjanes_columnes(observacions):
    totes_les_files = []
    for obs in observacions:
        totes_les_files.extend(obs)
    
    matriu = np.array(totes_les_files)
    mitjanes = np.mean(matriu, axis=0)
    return mitjanes

columnes = ["x(0)", "φ", "λ", "μ", "σ", "δ", "fitness", "distància", "temps"]
execucions_totals=50
observacions_globals_SA = llegir_observacions("Execucions_SA_globals.txt", "Execucions SA")
mitjanes_globals_SA = calcular_mitjanes_columnes(observacions_globals_SA)

observacions_globals_GA = llegir_observacions("Execucions_GA_globals.txt","Execucions GA continu")
mitjanes_globals_GA = calcular_mitjanes_columnes(observacions_globals_GA)

print("Proporció d'èxit SA: {}\n".format(len(observacions_globals_SA)/execucions_totals))
print("Mitjana dels òptims globals del SA:\n")

for i, mitjana in enumerate(mitjanes_globals_SA):
    print(f"{columnes[i]}: {mitjana:.10f}")

print("\n\n")
print("Proporció d'èxit GA: {}\n".format(len(observacions_globals_GA)/execucions_totals))
print("Mitjana dels òptims globals del SA:\n")

for i, mitjana in enumerate(mitjanes_globals_GA):
    print(f"{columnes[i]}: {mitjana:.10f}")
