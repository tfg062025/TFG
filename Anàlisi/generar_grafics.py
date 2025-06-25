import os
import re
import matplotlib.pyplot as plt

carpeta_fitxers = os.path.join(os.path.dirname(__file__), '..', 'Resultats', 'Fitxers evolució fitness GA continu')
carpeta_sortida = os.path.join(os.path.dirname(__file__), '..', 'Figures', 'Evolució del fitness millor individu GA continu')
os.makedirs(carpeta_sortida, exist_ok=True)

dades_per_T = {15: {}, 75: {}, 125: {}}
fitness_star=2566.999667640135158
nom = r'GA_continu_evolucio_P15000_PM([0-9.]+)_T(15|75|125)'

for nom_fitxer in os.listdir(carpeta_fitxers):
    hiperparam = re.match(nom, nom_fitxer)
    if hiperparam:
        pm = float(hiperparam.group(1))
        T = int(hiperparam.group(2))
        ruta_fitxer = os.path.join(carpeta_fitxers, nom_fitxer)
        with open(ruta_fitxer, 'r') as f:
            valors = [float(fila.strip()) for fila in f if fila.strip()]
            dades_per_T[T][pm] = valors


def estadistiques(T):
    print("T={}\n".format(T))
    for pm, valors in sorted(dades_per_T[T].items()):
        print("P_M={} \n".format(pm))
        for i, fit in enumerate(dades_per_T[T][pm]):
            if abs(fit - fitness_star) < 0.1:
                gen = i
                break

                    
        print("Generació d'aproximació a l'optim:{}".format(gen))
        print("Generacions necessàries per convergir:{}\n".format(len(dades_per_T[T][pm])))
    
def dibuixar_per_T(T, dades):
    
    for pm, valors in sorted(dades.items()):
        generacions = list(range(len(valors)))
        linia = plt.plot(generacions, valors, label=f'$P_M$ = {pm}')
        color = linia[0].get_color()
        plt.plot(generacions[0], valors[0], 'o', color=color, zorder=10)  
        plt.plot(generacions[-1], valors[-1], 's', color=color, zorder=10)  
        plt.legend(fontsize=50)  
        plt.xlabel('Generació')
        plt.ylabel('Fitness')
        plt.legend()


    nom_sortida = os.path.join(carpeta_sortida, f'grafica_T{T}.png')
    plt.savefig(nom_sortida, dpi=300)
    plt.show()

for T in [15, 75, 125]:
    dibuixar_per_T(T, dades_per_T[T])
    estadistiques(T)
