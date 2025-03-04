import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np
from Read_files import *

# from numba import njit
# from DBconn import *
from DatabaseConn import *


def seuil_critere(s_evalue, s_id, s_cover, b_evalue, b_id, b_cover):
    # s pour seuil, b pour blast
    condition = True
    if b_evalue > s_evalue or b_id < s_id or b_cover < s_cover:
        condition = False
    return condition


##########################################################################
## FONCTIONS CONDITIONS
def criteres(var_crit, seuil_crit):
    """ Retourne True si le gene passe le critère, sinon False
    var_crit : 0 = identite, 1 = evalue, 2 = couverture"""
    bool_crit = False
    match var_crit:  # 0 : identite, 1 : evalue, 2, couverture
        case 0:
            print(seuil_crit, "ident")
            bool_crit = True
        case 1:
            print(seuil_crit, "eval")
            bool_crit = True
        case 2:
            print(seuil_crit, "cover")
            bool_crit = True
    return bool_crit


def matrice_dotplot(list_hits, len_genome1, len_genome2):
    dotplot = np.zeros((len_genome1, len_genome2))
    print(len_genome1, len_genome2)
    for i in range(len(list_hits)):
        pos1, pos2 = list_hits[i][2], list_hits[i][3]
        dotplot[pos1-1][pos2-1] = 1
    return dotplot


# @njit
def threshold(mat_window,fenetre, stringence):
    count   = 0
    count2  = 0
    valid   = False
    for i in range(fenetre):
        if mat_window[i][i] == 1:
            count +=1
        if mat_window[i][fenetre-i-1]:
            count2 +=1
    if count >= stringence or count2 >= stringence:
        valid = True
    return valid

# @njit
def doplot_fenetre(matrice, fenetre, stringence):
    """Retrourne deux listes contenant les positions (coords x et y) des hits passant la stringence"""
    print("ENTRE")
    list1, list2 = [], []
    for i in range(len(matrice) - fenetre):
        for j in range(len(matrice[0]) - fenetre):
            window = matrice[i:i+fenetre, j:j+fenetre]
            if threshold(window, fenetre, stringence):
                list1.append(i)
                list2.append(j)
    return list1, list2


def dotplot_final(name_blast, type_seuil, val_seuil, val_fenetre, val_stringence):
    list_hits, len_query, len_subject = get_hits(name_blast, type_seuil, val_seuil)
    matrix_hits     = matrice_dotplot(list_hits, len_query, len_subject)
    list1, list2 =  doplot_fenetre(matrix_hits, val_fenetre, val_stringence)
    return list1, list2, len_query, len_subject
