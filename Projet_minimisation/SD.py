__authors__ = ("Ragousandirane Radjasandirane", "Kevin Mosca")
__version__ = "1.0.0"
__date__ = "2021-03-14"
__license__ = "BSD 3-Clause"

import math
import copy
import matplotlib.pyplot as plt
import sys

N = 6  # Nombre d'atomes total

K_BOND = 450
K_THETA = 55

THETA_EQ = 1.82421813  # 104.5200 EN RADIAN
L_EQ = 0.9572
R_MIN_HH = 0.449
R_MIN_OO = 3.5364
R_MIN_OH = 1.9927
EPSILON_HH = 0.046
EPSILON_OO = 0.1521
EPSILON_OH = 0.0836
Q_OH = -0.834
Q_H = 0.417
EPSILON_EQ = 1

# Indices
OH2 = 0
H1 = 1
H2 = 2
OH2_2 = 3
H1_2 = 4
H2_2 = 5

indice = {"OH2" :0 ,"H1" :1 ,"H2" :2 ,"OH2_2" : 3 ,"H1_2" :4 ,"H2_2" :5}

def evaluate_pdb(file_name):
    """
    test la validité d'un fichier .pdb pour une minimisation.

    Parameters
    ----------
    filename : str
        Nom du fichier pdb à tester.

    Returns
    -------
    True : boolean
        Booléen renvoyé pour valider le fichier d'entrée.
    """
    try:
        open(file_name)
    except ValueError:
        print("Le fichier de données n'existe pas")
        exit()
    if file_name[-4:] != ".pdb":
        print("Pas le bon format de données en entrée.")
        exit()
    with open(file_name, "r") as filin:
        lines = filin.readlines()
        for line in lines:
            if line.startswith("ATOM"):
                items = line.split()
                if len(items) != 10:
                    print("Pas le bon format de données en entrée.")
                    exit()
                try:
                    int(items[4])
                except ValueError:
                    print("Pas le bon format de données en entrée.")
                    exit()
                try:
                    float(items[5])
                except ValueError:
                    print("Pas le bon format de données en entrée.")
                    exit()
                try:
                    float(items[6])
                except ValueError:
                    print("Pas le bon format de données en entrée.")
                    exit()
                try:
                    float(items[7])
                except ValueError:
                    print("Pas le bon format de données en entrée.")
                    exit()

    return True

def read_pdb(file_name):
    """
    Lis les coordonnées dans un fichier .pdb.

    Parameters
    ----------
    filename : str
        Nom du fichier pdb.

    Returns
    -------
    atoms : list
        Liste de dictionnaire d'atomes.
    """
    atoms = []
    with open(file_name, "r") as filin:
        lines = filin.readlines()
        for line in lines:
            if line.startswith("ATOM"):
                items = line.split()
                atoms.append({"name": items[2],
                              "mol": items[4],
                              "x": float(items[5]),
                              "y": float(items[6]),
                              "z": float(items[7])})
    return atoms


def write_pdb(atoms, file_name):
    """
    Écrit les coordonnées minimisées dans un fichier.pdb.

    Parameters
    ----------
    atoms : list
        Liste de dictionnaire d'atomes.
    filename : str
        Nom du fichier pdb à créer.
    """
    i = 0
    with open(file_name, "w") as filin:
        for dico in atoms:
            i += 1
            coords = "{:8.3f}{:8.3f}{:8.3f}".format(dico["x"],
                                                    dico["y"],
                                                    dico["z"])
            filin.write("{:6s}{:5d} {:<3s}  {:3s} {:>4s}    {}{:6.2f}{:6.2f}\n"
            .format("ATOM", i, dico["name"],
            "TIP3", dico["mol"], coords, 1.00, 0.00))

    print("\nFile {} writed !\n".format(file_name))

def calc_en(atoms, Print=0):
    """
    Calcule de l'énergie potentielle.

    Parameters
    ----------
    atoms : list
        Liste de dictionnaire d'atomes.
    Print : int
        Indicatif pour écriture ou non des informations d'énergies

    Returns
    -------
    sum(Epot) : float
        Les énergies potentielles de chaque molécules sommées
    """
    Epot = []
    Nmol = int(atoms[-1]["mol"])

    # Calcul des distances

    # Entre OH2 et H1
    l_h1 = math.sqrt((atoms[OH2]["x"] - atoms[H1]["x"])**2 +
                     (atoms[OH2]["y"] - atoms[H1]["y"])**2 +
                     (atoms[OH2]["z"] - atoms[H1]["z"])**2)

    # Entre OH2 et H2
    l_h2 = math.sqrt((atoms[OH2]["x"] - atoms[H2]["x"])**2 +
                     (atoms[OH2]["y"] - atoms[H2]["y"])**2 +
                     (atoms[OH2]["z"] - atoms[H2]["z"])**2)

    # Énergie de liaison
    el_h1 = K_BOND * (l_h1 - L_EQ)**2
    el_h2 = K_BOND * (l_h2 - L_EQ)**2

    # Coordonnées des vecteurs pour le calcul d'angle
    X_H1 = atoms[H1]["x"] - atoms[OH2]["x"]
    Y_H1 = atoms[H1]["y"] - atoms[OH2]["y"]
    Z_H1 = atoms[H1]["z"] - atoms[OH2]["z"]

    X_H2 = atoms[H2]["x"] - atoms[OH2]["x"]
    Y_H2 = atoms[H2]["y"] - atoms[OH2]["y"]
    Z_H2 = atoms[H2]["z"] - atoms[OH2]["z"]

    # Énergie d'angle
    angle_h1_h2 = math.acos((X_H1 * X_H2 + Y_H1 * Y_H2 +
                             Z_H1 * Z_H2) /
                            (l_h1 * l_h2))
    theta_h1_2 = K_THETA * (angle_h1_h2 - THETA_EQ)**2

    Epot.append(el_h1 + el_h2 + theta_h1_2)

    # Si on a deux molécules, on calcul les autres termes
    if Nmol == 2:
        l_h1_2 = math.sqrt((atoms[OH2_2]["x"] - atoms[H1_2]["x"])**2 +
                           (atoms[OH2_2]["y"] - atoms[H1_2]["y"])**2 +
                           (atoms[OH2_2]["z"] - atoms[H1_2]["z"])**2)

        l_h2_2 = math.sqrt((atoms[OH2_2]["x"] - atoms[H2_2]["x"])**2 +
                           (atoms[OH2_2]["y"] - atoms[H2_2]["y"])**2 +
                           (atoms[OH2_2]["z"] - atoms[H2_2]["z"])**2)

        el_h1_2 = K_BOND * (l_h1_2 - L_EQ)**2
        el_h2_2 = K_BOND * (l_h2_2 - L_EQ)**2

        X_H1_2 = atoms[H1_2]["x"] - atoms[OH2_2]["x"]
        Y_H1_2 = atoms[H1_2]["y"] - atoms[OH2_2]["y"]
        Z_H1_2 = atoms[H1_2]["z"] - atoms[OH2_2]["z"]

        X_H2_2 = atoms[H2_2]["x"] - atoms[OH2_2]["x"]
        Y_H2_2 = atoms[H2_2]["y"] - atoms[OH2_2]["y"]
        Z_H2_2 = atoms[H2_2]["z"] - atoms[OH2_2]["z"]
        angle_h1_h2_2 = math.acos((X_H1_2 * X_H2_2 +
                                              Y_H1_2 * Y_H2_2 +
                                              Z_H1_2 * Z_H2_2) /
                                             (l_h1_2 * l_h2_2))
        theta_h1_2_2 = K_THETA * (angle_h1_h2_2 - THETA_EQ)**2

        # Distances entre atomes
        dist_inter = []
        list_atomes = ["H1_H1_2","H1_OH2_2","H1_H2_2","OH2_H1_2","OH2_OH2_2","OH2_H2_2","H2_H1_2","H2_OH2_2","H2_H2_2"]
        for l in list_atomes:
            items = l.split("_", maxsplit = 1)
            if items[0][0] == "H":
                q1 = Q_H
                if items[1][0] == "H":
                    epsi = EPSILON_HH
                    rmin = R_MIN_HH
                    q2 = Q_H
                elif items[1][0] == "O":
                    epsi = EPSILON_OH
                    rmin = R_MIN_OH
                    q2 = Q_OH
            else:
                q1 = Q_OH
                if items[1][0] == "H":
                    epsi = EPSILON_OH
                    rmin = R_MIN_OH
                    q2 = Q_H
                elif items[1][0] == "O":
                    epsi = EPSILON_OO
                    rmin = R_MIN_OO
                    q2 = Q_OH
            dic = {"r": math.sqrt((atoms[indice[items[0]]]["x"] - atoms[indice[items[1]]]["x"])**2 +
                                  (atoms[indice[items[0]]]["y"] - atoms[indice[items[1]]]["y"])**2 +
                                  (atoms[indice[items[0]]]["z"] - atoms[indice[items[1]]]["z"])**2),
                   "epsilon": epsi, "r_min": rmin,
                   "q1": q1, "q2": q2}
            dist_inter.append(dic)

        # Énergie de Lennard-Jones et Coulomb
        coulomb = 0
        lennard_jones = 0
        for i in dist_inter:
            lennard_jones += i["epsilon"]*((i["r_min"]/i["r"])**12 -
                                           2*(i["r_min"]/i["r"])**6)
            coulomb += ((i["q1"]*i["q2"])/( i["r"])) * 332.0716
        Epot.append(lennard_jones + coulomb + el_h1_2 +
                    el_h2_2 + theta_h1_2_2)

    # Affiche les valeurs pour une molécule
    if Print:
        print("Valeur pour la molécule 1 : ")
        print("Liaison OH2-H1 : {}\nLiaison OH2-H2 : {}".format(l_h1, l_h2))
        print("Angle : {}\n".format(angle_h1_h2))

        if Nmol == 2:   
            print("Valeur pour la molécule 2 : ")
            print("Liaison OH2-H1 : {}\nLiaison OH2-H2 : {}".format(l_h1_2, l_h2_2))
            print("Angle : {}".format(angle_h1_h2_2))
            print("LJ : {}".format(lennard_jones))
            print("Coulomb : {}\n".format(coulomb))
            


    return sum(Epot)


def derivate(atoms, delta):
    """
    Dérivation numérique de l'énergie selon les 9 variables.

    Parameters
    ----------
    atoms : list
        Liste de dictionnaire d'atomes.
    delta : float
        Valeur pour la dérivé numérique

    Returns
    -------
    grad : lst
        Liste de dictionnaire contenant le gradient des dérivés numériques.
    """
    grad = []

    for i in range(len(atoms)):
        # Copie des positions de toutes les molécules
        # pour effectuer les variations de delta
        atoms_x_pdelta = copy.deepcopy(atoms)
        atoms_y_pdelta = copy.deepcopy(atoms)
        atoms_z_pdelta = copy.deepcopy(atoms)
        atoms_x_mdelta = copy.deepcopy(atoms)
        atoms_y_mdelta = copy.deepcopy(atoms)
        atoms_z_mdelta = copy.deepcopy(atoms)

        atoms_x_pdelta[i]["x"] += delta
        atoms_y_pdelta[i]["y"] += delta
        atoms_z_pdelta[i]["z"] += delta
        atoms_x_mdelta[i]["x"] -= delta
        atoms_y_mdelta[i]["y"] -= delta
        atoms_z_mdelta[i]["z"] -= delta

        # Gradient sous forme de liste de dict
        x = (calc_en(atoms_x_pdelta) - calc_en(atoms_x_mdelta))/(2*delta)
        y = (calc_en(atoms_y_pdelta) - calc_en(atoms_y_mdelta))/(2*delta)
        z = (calc_en(atoms_z_pdelta) - calc_en(atoms_z_mdelta))/(2*delta)
        grad.append({"x": x,
                     "y": y,
                     "z": z})
    return grad


def calc_GRMS(grad):
    """
    Calcule le GRMS à chaque étape de la minimisation.

    Parameters
    ----------
    atoms : list
        Liste de dictionnaire d'atomes.
    delta : float
        Valeur pour la dérivé numérique

    Returns
    -------
    math.sqrt((som)/(N*3)) : float
        Valeur du GRMS
    """
    som = 0

    # Somme de toutes les valeurs du gradient
    # en parcourant les dico de la liste
    for dico in grad:
        som += dico["x"]**2
        som += dico["y"]**2
        som += dico["z"]**2

    return math.sqrt((som)/(N*3))


def calc_decent(atoms):
    """
    Calcule du pas prochain selon l'opposé du gradient tant que le GRMS et < à 0.01.

    Parameters
    ----------
    atoms : list
        Liste de dictionnaire d'atomes.
    """
    threshold = 0.01
    Lambda = 1e-04
    DELTA = Lambda / 10

    # Premier calcul du gradient au pas inital
    step = 1
    grad = derivate(atoms, DELTA)
    GRMS = calc_GRMS(grad)
    energie = calc_en(atoms)

    # Initalisation des listes avec la 1ere valeur
    energie_value = [energie]
    GRMS_value = [GRMS]

    print("Energie : ", calc_en(atoms))
    print("GRMS : ", GRMS)

    print("\n")
    while (GRMS > threshold) and step != 10000:
        step += 1

        # Un pas de gradient
        for i in range(len(atoms)):
            atoms[i]["x"] -= Lambda * grad[i]["x"]
            atoms[i]["y"] -= Lambda * grad[i]["y"]
            atoms[i]["z"] -= Lambda * grad[i]["z"]
        # Calcul du gradient avec les nouvelles positions
        grad = derivate(atoms, DELTA)
        GRMS = calc_GRMS(grad)
        energie = calc_en(atoms)
        energie_value.append(energie)
        GRMS_value.append(GRMS)

        print("Energie : ", energie)
        print("GRMS :", GRMS)
        print("STEP : ", step)
        print("\n")

    calc_en(atoms, Print=1)
    x = step/2
    y = max(max(GRMS_value), max(energie_value)) / 1.7
    step_value = list(range(step))

    plt.plot(step_value, energie_value, label="Energie")
    plt.plot(step_value, GRMS_value, label="GRMS")

    plt.title("Energie et GRMS en fonction de step")
    plt.xlabel("STEP")
    plt.ylabel("Energie et GRMS")
    plt.text(x, y, "GRMS : {} \nEnergie : {}".format(GRMS, energie))
    plt.legend(loc="best")

    # Ecriture de deux fichier (changer le nom selon les coords utilisées)
    write_pdb(atoms, sys.argv[2])
    plt.savefig("graph_{}.png".format(sys.argv[2].split(".")[0]))
    plt.show()

    print(" \n\nNOMBRE DE PAS : ", step)


if __name__ == "__main__":
    # Coords du TP MINI avec une seule molécule d'eau
    atoms_TP_mini = [{"name": "OH2", "mol": "1",
                      "x": 0.931, "y": 0.293, "z": 1.371},
                     {"name": "H1", "mol": "1",
                      "x": 1.568, "y": 1.002, "z": 1.273},
                     {"name": "H2", "mol": "1",
                      "x": 1.464, "y": -0.488, "z": 1.523}]

    # Coords du TP 1
    atoms_TP_1 = [{"name": "OH2", "mol": "1",
                   "x": -0.765, "y": -0.097, "z": 0.01},
                  {"name": "H1", "mol": "1",
                   "x": -0.222, "y": 0.691, "z": 0.024},
                  {"name": "H2", "mol": "1",
                   "x": -0.138, "y": -0.819, "z": -0.035}]

    # Coords d'une molécule d'eau plane
    atoms_plan = [{"name": "OH2", "mol": "1",
                   "x": -0.765, "y": -0.097, "z": 0.01},
                  {"name": "H1", "mol": "1",
                   "x": -1.765, "y": -0.097, "z": 0.01},
                  {"name": "H2", "mol": "1",
                   "x": 0.235, "y": -0.097, "z": 0.01}]

    # Coords de deux molécules d'eau
    atoms_deux_eaux = [{"name": "OH2", "mol": "1",
                        "x": 0.931, "y": 0.293, "z": 1.371},
                       {"name": "H1", "mol": "1",
                        "x": 1.568, "y": 1.002, "z": 1.273},
                       {"name": "H2", "mol": "1",
                        "x": 1.464, "y": -0.488, "z": 1.523},
                       {"name": "OH2", "mol": "2",
                        "x": -1.279, "y": 0.298, "z": 1.358},
                       {"name": "H1", "mol": "2",
                        "x": -1.860, "y": 1.049, "z": 1.476},
                       {"name": "H2", "mol": "2",
                        "x": -1.820, "y": -0.461, "z": 1.575}]

    atoms_deux_eaux_TP1 = [{"name": "OH2", "mol": "1",
                            "x": -0.765, "y": -0.097, "z": 0.01},
                           {"name": "H1", "mol": "1",
                            "x": -0.222, "y": 0.691, "z": 0.024},
                           {"name": "H2", "mol": "1",
                            "x": -0.138, "y": -0.819, "z": -0.035},
                           {"name": "OH2", "mol": "2",
                            "x": -3.185, "y": -0.039, "z": -0.928},
                           {"name": "H1", "mol": "2",
                            "x": -3.544, "y": 0.829, "z": -1.113},
                           {"name": "H2", "mol": "2",
                            "x": -3.921, "y": -0.637, "z": -1.058}]
    if len(sys.argv) != 3:
        print("Mauvais appel du script")
        exit()
    elif evaluate_pdb(sys.argv[1]) and sys.argv[2][-4:] == ".pdb":
        calc_decent(read_pdb(sys.argv[1]))

