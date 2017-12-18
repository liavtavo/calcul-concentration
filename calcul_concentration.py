
# coding: utf-8

# In[7]:



# Programme pour calculer la masse moléculaire à partir d'une formule brute,
# et calculer la composition d'une solution de plusieurs composants.
# Programme élaboré dans le cadre du DU Bioinformatique et Datamanagement
# de Ioannis Nicolis, Université Paris Descartes. Année 2017-2018.
# Auteur : Thomas Liautaud

# définition d'une liste d'exemples :
exemples = ['NaCl',
            'Fe2(SO4)3.7H2O',
            'K4[Fe(CN)6]',
            'Ag3AsS3',
            'C12H19NO2',
            'Ca2NaH(SiO3)3',
            'Fe2Al9O6(SiO4)4(OOH)2',
            'CH3[CH2]6CH3',
            'CH3CH(CH3)CH2CH3',
            '[FeCl2(H2O)4]Cl.2H2O',
            'C3H3FeO6',
            'Cu5(SiO3)4(OH)2',
            'CH3[CH2]12CH3',
            'Cr3Fe2O12'
           ]

# définition d'une liste d'indices :
l = list(range(0,201))
indices = list()
for i in l:
    indices.append(str(i))
    
# définition des séparateurs nécessaires à la décomposition des formules brutes :
majuscules = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
minuscules = majuscules.lower()
parentheses  = '()'
crochets = '[]'
delimiteurs = [majuscules, parentheses, crochets, '.',indices]
non_del = [minuscules]

# liste des masses molaires des atomes :
mass = {'H': 1.00794, 'He': 4.002602, 'Li': 6.941, 'Be': 9.012182,
        'B': 10.811, 'C': 12.0107, 'N': 14.0067, 'O': 15.9994,
        'F': 18.9984032, 'Ne': 20.1797, 'Na': 22.98976928,
        'Mg': 24.305, 'Al': 26.9815386, 'Si': 28.0855,
        'P': 30.973762, 'S': 32.065, 'Cl': 35.453, 'Ar': 39.948,
        'K': 39.0983, 'Ca': 40.078, 'Sc': 44.955912,
        'Ti': 47.867, 'V': 50.9415, 'Cr': 51.9961, 'Mn': 54.938045,
        'Fe': 55.845, 'Co': 58.933195, 'Ni': 58.6934, 'Cu': 63.546,
        'Zn': 65.409, 'Ga': 69.723, 'Ge': 72.64, 'As': 74.9216,
        'Se': 78.96, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.4678,
        'Sr': 87.62, 'Y': 88.90585, 'Zr': 91.224, 'Nb': 92.90638,
        'Mo': 95.94, 'Tc': 98.0, 'Ru': 101.07, 'Rh': 102.9055,
        'Pd': 106.42, 'Ag': 107.8682, 'Cd': 112.411, 'In': 114.818,
        'Sn': 118.71, 'Sb': 121.76, 'Te': 127.6, 'I': 126.90447,
        'Xe': 131.293, 'Cs': 132.9054519, 'Ba': 137.327, 'La': 138.90547,
        'Ce': 140.116, 'Pr': 140.90765, 'Nd': 144.242, 'Pm': 145.0,
        'Sm': 150.36, 'Eu': 151.964, 'Gd': 157.25, 'Tb': 158.92534,
        'Dy': 160.5, 'Ho': 164.93032, 'Er': 167.259, 'Tm': 168.93421,
        'Yb': 173.04, 'Lu': 174.967, 'Hf': 178.49, 'Ta': 180.94788,
        'W': 183.84, 'Re': 186.207, 'Os': 190.23, 'Ir': 192.217,
        'Pt': 195.084, 'Au': 196.966569, 'Hg': 200.59, 'Tl': 204.3833,
        'Pb': 207.2, 'Bi': 208.9804, 'Po': 209.0, 'At': 210.0,
        'Rn': 222.0, 'Fr': 223.0, 'Ra': 226.0254, 'Ac': 227.0,
        'Th': 232.03806, 'Pa': 231.03588, 'U': 238.02891, 'Np': 237.0,
        'Pu': 244.06, 'Am': 241.06, 'Cm': 247.0, 'Bk': 247.0, 'Cf': 251.0,
        'Es': 252.0, 'Fm': 257.0, 'Md': 258.0, 'No': 259.0, 'Lr': 262.0,
        'Rf': 261.0, 'Db': 262.0, 'Sg': 266.0, 'Bh': 264.0,
        'Hs': 269.0, 'Mt': 268.0, 'Ds': 281.0, 'Rg': 280.0,
        'Cn': 285.0, 'Fl': 289.0
       }

# Définition d'une fonction permettant de concatener les chiffres consécutifs dans une liste :
def concatenation_ind(ch):
    """
    Fonction permettant de concatener les chiffres consécutifs dans une liste.
    """
    chiffres = '0123456789'
    i = 1
    while i < len(ch)-1:
        if ch[i] in chiffres and ch[i+1] in chiffres:
            ch[i:i+2] = [ch[i]+ch[i+1]]
        i+=1
    return ch

# Définition d'une fonction permettant de de retourner un dictionnaire des éléments et de leurs nombres,
# i.e. une formule brute sous forme de dictionnaire:
def denombrement(decompo):
    """
    fonction permettant de de retourner un dictionnaire des éléments
    et de leurs nombres, i.e. une formule brute sous forme de dictionnaire.
    """
    nombre = {}
    l = decompo.split()
    # utilise la fonction définie précédemment :
    concatenation_ind(l)
    while len(l) > 1:
        e, n, *f = l
        if n in indices:
            if e in nombre.keys():
                nombre[e] += float(n)
            else:
                nombre[e] = float(n)
            l = f
        else:
            if e in nombre.keys():
                nombre[e] += 1
            else:
                nombre[e] = 1
            l = [n] + f
    else:
        for e in l:
            if e in nombre.keys():
                nombre[e] += 1
            else:
                nombre[e] = 1
    return nombre

# Cette fonction permet à partir d'un dictionnaire de retourner une chaine de caractère de type str,
# en duplicant les clés du dictionnaire autant de fois que la valeur associée à la clée :
def concatenation(dico):
    """
    Cette fonction permet à partir d'un dictionnaire de retourner une chaine de caractère de type str,
    en duplicant les clés du dictionnaire autant de fois que la valeur associée à la clée."
    """
    formule = ''
    for k in dico.keys():
        formule = formule + k*int(dico[k])
    return formule

# fonction pour encadrer les caractères entre parenthèses par des espaces
# et les indices des groupes par des espaces :
def enlev_parenthese(formule):
    """
    fonction pour encadrer les caractères entre parenthèses par des espaces
    et les indices des groupes par des espaces.
    NB : ne fonctionne pas si indices > 9
    """
    f = formule
    g = ''
    while len(f) > 0:
        if f[0] in ['(',')']:
            g = g + ' '
            try:
                if f[1] in indices:
                    g = g + f[1] + ' '
                    f = f[2:]
                else:
                    f = f[1:]
            except:
                f = f[1:]
        else:
            g = g + f[0]
            f = f[1:]
    return g

# Définition d'une fonction permettant de dupliquer les groupes entre crochets
# dans une formule brute autant de fois que l'indique l'indice du groupe :
def enlev_crochet(formule):
    """
    Fonction permettant de dupliquer les groupes entre crochets
    dans une formule brute autant de fois que l'indique l'indice du groupe.
    NB : ne fonctionne pas si indices > 9
    """
    f = formule
    g = ''
    while len(f) > 0:
        if f[0] in ['[',']']:
            g = g + ' '
            try:
                if f[1] in indices:
                    g = g + f[1] + ' '
                    f = f[2:]
                else:
                    f = f[1:]
            except:
                f = f[1:]
        else:
            g = g + f[0]
            f = f[1:]
    return g

# Définition d'une fonction permettant de dupliquer les composés d'addition (hydrates par ex.)
# au bout d'une chaine de type str, autant de fois que leur indice l'indique :
def hydrate(formule):
    """
    Définition d'une fonction permettant de dupliquer les composés d'addition (hydrates par ex.)
    au bout d'une chaine de type str, autant de fois que leur indice l'indique.
    NB : ne fonctionne pas si il y a deux composés d'addition consécutifs.
    """
    f = formule
    g = ''
    while len(f) > 0:
        if f[0] == '.':
            fin, f = f[1:], ''            
            facteur = ''
            while fin[0] in '0123456789':
                    facteur = facteur + fin[0]
                    fin = fin[1:]
            if facteur == '':
                g = g + fin
            else:
                g = g + fin*int(facteur)
        else:
            g = g + f[0]
            f = f[1:]
    return g

# Cette fonction permet de décomposer une formule brute de type str dans une liste,
# en séparant les atomes et les indices :
def decompo(formule):
    """
    Cette fonction permet de décomposer une formule brute de type str dans une liste,
    en séparant les atomes et les indices.
    """
    decompo = str()
    for element in formule:
        for sep in delimiteurs:
            if element in sep:
                decompo += ' ' + element
        for sepbis in non_del:
            if element in sepbis:
                decompo += element
    return decompo

# Cette fonction permet de calculer la masse molaire d'une molécule
# à partir d'une formule brute sous forme de dictionnaire :
def massmol(formule):
    """
    Cette fonction permet de calculer la masse molaire d'une molécule,
    à partir d'une formule brute sous forme de dictionnaire.
    """
    mm = 0.
    for atome in formule.keys():
        mm = mm + formule[atome]*mass[atome]
    return mm

# Maintenant que les fonctions sont définies, les lignes ci-dessous donnent :
# les invites pour saisir les composants, leurs formules et les quantités à traiter;
# le processus de traitement des données ;
# le traitement des erreurs;
# l'affichage des résultats :
composants = {}
invite = 'o'
while invite != 'n':
    if invite not in ['o', 'n']:
        print(f"Je n'ai pas compris la réponse : {invite}")
        print('---------------')
        invite = input("Voulez-vous saisir un nouveau composant (o/n) ?")
    else:
        nom = input("Saisissez le nom du composant : ")
        f = input("Saisissez la formule brute ou semi-développée : ")
        try:
            r = enlev_parenthese(f)
            s = denombrement(r)
            t = concatenation(s)
            u = enlev_crochet(t)
            v = denombrement(u)
            w = concatenation(v)
            x = hydrate(w)
            y = decompo(x)
            z = denombrement(y)
            a = massmol(z)
            composants[nom] = {'formule': f}
            composants[nom]['massmol'] = a
            print(f"""Vous avez saisi le composant "{nom}" de formule brute "{composants[nom]['formule']}" """)
            confirmation = input("Confirmez-vous la saisie (o/n) ? ")
            while confirmation not in ['o', 'n']:
                print(f"""Je n'ai pas compris la réponse "{confirmation}" """)
                print("veuillez répondre par o ou n.")
                confirmation = input("Confirmez-vous la saisie (o/n) ? ")
            if confirmation == 'n':
                del composants[nom]
                invite = 'o'         
            else:
                print(f"La masse molaire de {nom} = {composants[nom]['massmol']:.3f} g/mol")
                quantite = -1
                while quantite < 0:
                    try:
                        print('---------------')
                        quantite = float(input("Saisissez la quantité exprimée en mg : "))
                        if quantite > 0:
                            composants[nom]['quantite'] = quantite
                            print(f"Vous avez saisi {composants[nom]['quantite']:.2f} mg de {nom}")
                            confirmation_q = input("Confirmez-vous la saisie (o/n) ? ")
                            while confirmation_q not in ['o', 'n']:
                                print(f"""Je n'ai pas compris la réponse "{confirmation_q}" """)
                                print("veuillez répondre par o ou n.")
                                confirmation_q = input("Confirmez-vous la saisie (o/n) ? ")
                            if confirmation_q == 'n':
                                del composants[nom]['quantite']
                                quantite = -1         
                            else:
                                print('---------------')
                                invite = input("Voulez-vous saisir un nouveau composant (o/n) ? ")
                        else:
                            print("veuillez saisir une valeur positive.")
                            quantite = -1
                    except ValueError:
                        print("Veuillez saisir une valeur numérique.")
                        quantite = -1
        except KeyError:
            print(f"""ATTENTION : La formule "{f}" est incorrecte.""")
            print(">>>> Utilisez exclusivement les caractères alphanumériques et les caractères suivants : .()[] ")
            print(">>>> Respecter la casse : 'NaCl' et non 'nacl' ou NACL")
            print("Exemples : ")
            for c in exemples:
                print(c)
            invite = 'o'

# Les lignes ci-dessous permettent de :
# saisir le volume d'eau;
# calculer les résultats demandés ;
# afficher les résultats :
confirmation_v = 'n'
while confirmation_v is not 'o':
    v = -1
    while v < 0:
        try:
            print('---------------')
            v = float(input("saisissez le volume d'eau en mL : "))
            if v > 0: 
                composants['eau'] = {'formule': 'H2O', 'massmol': 18, 'quantite': 1000*v}
                print(f"vous avez saisi {v} mL d'eau")
                confirmation_v = input("confirmez vous le volume saisi (o/n) : ")
                if confirmation_v not in ['o', 'n']:
                    print("Veuillez répondre par o ou n.")
                    confirmation_v = input("confirmez vous le volume saisi (o/n) : ")
                elif confirmation_v is'n':
                    del composants['eau']['quantite']
            else:
                print("Veuillez saisir une valeur positive.")
                v = -1
        except ValueError:
            print("Veuillez saisir une valeur numérique.")
            v = -1

masse_composants = 0
for comp in composants.keys():
    masse_composants = masse_composants + composants[comp]['quantite']
for comp2 in composants.keys():
    Cmg = 1000*composants[comp2]['quantite']/composants['eau']['quantite']
    composants[comp2]['Cg/g'] = Cmg
    Cmol = ((composants[comp2]['quantite']*1000)/composants[comp2]['massmol'])*(1000/composants['eau']['quantite'])
    composants[comp2]['Cmol/L'] = Cmol
    fmass = 1000*composants[comp2]['quantite']/masse_composants
    composants[comp2]['fmass'] = fmass
for c in composants.keys():
    print('---------------')
    print(f"composant : {c}")
    print(f"formule brute : {composants[c]['formule']}")
    print(f"masse molaire : {composants[c]['massmol']:.3f} g/mol")
    print(f"fraction massique : {composants[c]['fmass']:.3f} mg de {c}/g de solution")
    print(f"concentration : {composants[c]['Cg/g']:.3f} g/L soit {composants[c]['Cmol/L']:.2f} mmol/L")
                    

