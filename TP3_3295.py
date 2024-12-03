#Laurent Bouchard 20184162 laurent.bouchard@umontreal.ca
#DaoHan Li 20209828 bloodthurster.sd@gmail.com

from Bio import SeqIO
import pandas as pd
import math

""" 
Section pour les fonctions
"""
#1-On fait une fonction pour lire les fichiers fasta
# retourne un df contenant toutes les informations du fichier
def read_fasta(file_path):
    data = []
    for record in SeqIO.parse(file_path, "fasta"):
        parts = record.description.split("|")
        id_part = parts[0]
        anticodon = parts[1]
        info_part = parts[2] if len(parts) > 2 else None
        seq = str(record.seq)
        data.append({"ID": id_part, "Anticodon": anticodon, "Info": info_part, "Sequence": seq})

    #On tranforme tout en df
    df = pd.DataFrame(data)
    return df


#2-On fait une fonction pour générer les kmers de taille k avec leur pos dans la seq
#On retourne la db avec un nouvelle col
"""Vieille Version
def generate_kmers(sequence, k):
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        kmers.append((kmer, i))
    return kmers
""" 
def generate_kmers_with_seed(sequence, seed):
    k = len(seed)
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        kmers.append((kmer, i))
    return kmers #, seed


#3-On recherche les kmers dans la base de donnée
# retourne les positions où il est trouvé
"""Vieille Version
def find_kmers_in_db(kmers, df_dataBase):
    #On initie la variable des match
    matches = set()
    #Pour tous les kmers d'un S, ...
    for kmer, query_pos in kmers:
        #On regarde toutes les rangées de la db
        for _, row in df_dataBase.iterrows():
            #On identifie la rangée
            db_id = row["ID"]
            db_anticodon = row["Anticodon"]
            #db_info = row["Info"]
            db_seq = row["Sequence"]
            #SI on trouve le kmer dans la S, ...
            if kmer in db_seq:
                #On trouve la position dans la db.
                db_pos = db_seq.index(kmer)
                matches.add((kmer, query_pos, db_id, db_anticodon, db_pos))
                
                #On semble faire plusieurs fois la même comparaison
                #Ou plutot, on compare les lettres
                #print("u kmer : " + kmer)
                #print("n kmer : " + db_seq[db_pos:db_pos+11])
                
    return list(matches)  # Convertir le set en liste
""" 
def find_kmers_with_seed(kmers, df_dataBase, seed):
    """
    Recherche des k-mers dans la base de données en respectant la graine.

    Args:
        kmers (list): Liste de tuples (k-mer, position).
        df_dataBase (DataFrame): Base de données contenant les séquences.
        seed (str): Graine (ex. "111010010100").

    Returns:
        list: Liste de correspondances sous forme (kmer, query_pos, db_id, db_anticodon, db_pos).
    """
    matches = []
    k = len(seed)
    for kmer, query_pos in kmers:
        for _, row in df_dataBase.iterrows():
            db_seq = row["Sequence"]
            db_id = row["ID"]
            db_anticodon = row["Anticodon"]

            # Rechercher des correspondances
            for i in range(len(db_seq) - k + 1):
                db_kmer = db_seq[i:i + k]
                # Vérifier la correspondance selon la graine
                match = all(
                    (kmer[j] == db_kmer[j]) if seed[j] == "1" else True
                    for j in range(k)
                )
                if match:
                    matches.append((kmer, query_pos, db_id, db_anticodon, i))
    return matches



#4- On extend HSP des deux côtés
def extend_hsp(hsp, query_seq, db_seq, E=4):
    #On récupère les données de hsp.
    query_pos, db_pos, kmer = hsp
    #Score de tous les nuclé
    max_score = len(kmer) * 5
    current_score = max_score

    aligned_query = kmer
    aligned_db = kmer

    # Étendre à gauche
    left_query_pos, left_db_pos = query_pos - 1, db_pos - 1
    while left_query_pos >= 0 and left_db_pos >= 0:
        # Vérifier que les indices sont dans les limites
        if left_query_pos < 0 or left_db_pos < 0:
            break

        q_char = query_seq[left_query_pos]
        d_char = db_seq[left_db_pos]
        aligned_query = q_char + aligned_query
        aligned_db = d_char + aligned_db

        # Calcul du score
        if q_char == d_char:
            current_score += 5  # Match
        else:
            current_score -= 4  # Mismatch

        # Mise à jour du score maximal
        max_score = max(max_score, current_score)

        # Arrêt si le score diminue trop
        if current_score < max_score - E:
            # Supprimer la dernière extension
            aligned_query = aligned_query[1:]
            aligned_db = aligned_db[1:]
            break

        left_query_pos -= 1
        left_db_pos -= 1

    # Réinitialiser les positions pour l'extension à droite
    right_query_pos, right_db_pos = query_pos + len(kmer), db_pos + len(kmer)
    current_score = max_score  # Repartir avec le score après l'extension gauche

    # Étendre à droite
    while right_query_pos < len(query_seq) and right_db_pos < len(db_seq):
        # Vérifier que les indices sont dans les limites
        if right_query_pos >= len(query_seq) or right_db_pos >= len(db_seq):
            break

        q_char = query_seq[right_query_pos]
        d_char = db_seq[right_db_pos]
        aligned_query += q_char
        aligned_db += d_char

        # Calcul du score
        if q_char == d_char:
            current_score += 5  # Match
        else:
            current_score -= 4  # Mismatch

        # Mise à jour du score maximal
        max_score = max(max_score, current_score)

        # Arrêt si le score diminue trop
        if current_score < max_score - E:
            # Supprimer la dernière extension
            aligned_query = aligned_query[:-1]
            aligned_db = aligned_db[:-1]
            break

        right_query_pos += 1
        right_db_pos += 1

    # Score final
    final_score = max_score
    return aligned_query, aligned_db, final_score


#5-Fonction pour appeler la fonction extend_hsp
def process_hsp_extensions(row, db_df, E=4):
    extensions = []
    #On traverse tous les matchs
    for match in row["Matches"]:
        kmer, query_pos, db_id, anticodon, db_pos = match
        
        #On récup la S de base
        db_seq = db_df.loc[
            (db_df["ID"] == db_id) & (db_df["Anticodon"] == anticodon),
            "Sequence"
        ].values[0]

        #On vérifie les indices
        if query_pos < 0 or query_pos >= len(row["Sequence"]) or \
           db_pos < 0 or db_pos >= len(db_seq):
            continue  # Ignorer les correspondances invalides
        
        # Étendre le HSP
        aligned_query, aligned_db, final_score = extend_hsp(
            (query_pos, db_pos, kmer),
            row["Sequence"],  # Séquence cible
            db_seq,           # Séquence de la base
            E
        )
        
        # Ajouter le résultat
        extensions.append((aligned_query, aligned_db, final_score))
    
    return extensions


#6- On fuse les HSP
def merge_overlapping_hsps(hsps):
    """
    Fusionne les HSPs chevauchants d'une paire de séquences.

    Args:
    - hsps (list of dict): Liste de HSPs sous forme de dictionnaires, chaque HSP contient :
        {
            "query_start": int,  # Position de début dans la séquence cible
            "query_end": int,    # Position de fin dans la séquence cible
            "db_start": int,     # Position de début dans la séquence de la base
            "db_end": int,       # Position de fin dans la séquence de la base
            "aligned_query": str, # Portion alignée de la séquence cible
            "aligned_db": str,    # Portion alignée de la séquence de la base
            "score": int          # Score de l'HSP
        }

    Returns:
    - list of dict: Liste des HSPs fusionnés.
    """
    if not hsps:
        return []

    # Trier les HSPs par position de début dans la séquence cible
    hsps = sorted(hsps, key=lambda h: h["query_start"])

    merged_hsps = [hsps[0]]

    for hsp in hsps[1:]:
        last = merged_hsps[-1]

        # Vérifier si les HSPs se chevauchent dans la séquence cible
        if hsp["query_start"] <= last["query_end"]:
            # Fusionner les intervalles
            new_query_start = min(last["query_start"], hsp["query_start"])
            new_query_end = max(last["query_end"], hsp["query_end"])
            new_db_start = min(last["db_start"], hsp["db_start"])
            new_db_end = max(last["db_end"], hsp["db_end"])

            # Fusionner les alignements
            """
            new_aligned_query = (
                last["aligned_query"][: hsp["query_start"] - last["query_start"]]
                + hsp["aligned_query"]
            )
            new_aligned_db = (
                last["aligned_db"][: hsp["db_start"] - last["db_start"]]
                + hsp["aligned_db"]
            )
            """
            new_aligned_query = last["aligned_query"][:hsp["query_start"] - last["query_start"]] + hsp["aligned_query"]
            new_aligned_db = last["aligned_db"][:hsp["db_start"] - last["db_start"]] + hsp["aligned_db"]
            overlap_len = last["query_end"] - hsp["query_start"] + 1
            if overlap_len > 0:
                new_aligned_query = (last["aligned_query"][: -overlap_len] + hsp["aligned_query"])
                new_aligned_db = (last["aligned_db"][: -overlap_len] + hsp["aligned_db"])

            # Combiner les scores
            new_score = last["score"] + hsp["score"]

            # Mettre à jour le dernier HSP
            merged_hsps[-1] = {
                "query_start": new_query_start,
                "query_end": new_query_end,
                "db_start": new_db_start,
                "db_end": new_db_end,
                "aligned_query": new_aligned_query,
                "aligned_db": new_aligned_db,
                "score": new_score,
            }
        else:
            # Ajouter un nouveau HSP si pas de chevauchement
            merged_hsps.append(hsp)

    return merged_hsps


#7- On fait une fonction pour convertir les extentions en dictionnaire
#Ca rend les choses plus clair.
def prepare_hsps(extensions):
    """
    Convertit les extensions en une liste de dictionnaires compatibles avec merge_overlapping_hsps.

    Args:
    - extensions (list of tuples): Liste des alignements étendus (aligned_query, aligned_db, final_score).

    Returns:
    - list of dict: Liste de HSPs sous forme de dictionnaires.
    """
    hsps = []
    for extension in extensions:
        aligned_query, aligned_db, final_score = extension

        hsp = {
            "query_start": 0,  # À ajuster si besoin d'utiliser des indices réels
            "query_end": len(aligned_query) - 1,
            "db_start": 0,  # À ajuster si besoin d'utiliser des indices réels
            "db_end": len(aligned_db) - 1,
            "aligned_query": aligned_query,
            "aligned_db": aligned_db,
            "score": final_score
        }
        hsps.append(hsp)
    return hsps


#8- On calc le bit score
"""Vieille Version
def calculate_bitscore_and_evalue(hsp, total_db_size, query_size, ss=1e-3):
    
    # Calcule le bitscore et la e-value d'un HSP.

    # Args:
    # - hsp (dict): Un HSP avec les champs "score" et autres.
    # - total_db_size (int): Taille totale des séquences dans la base de données.
    # - query_size (int): Taille de la séquence cible.
    # - ss (float): Seuil significatif pour la e-value.

    # Returns:
    # - dict: HSP enrichi avec les champs "bitscore", "evalue", et "significant".
    

    S = hsp["score"]  # Score brut
    K = 0.192  # Constante
    ln2 = math.log(2)  # ln(2)

    # Calcul du bitscore
    try:
        B = round((S * math.log(K)) / ln2)
    except ValueError:
        B = 0  # Si un score brut est invalide, on retourne 0 par défaut

    # Calcul de la e-value avec ajustements logarithmiques
    ln_total_db_size = math.log(total_db_size) if total_db_size > 0 else 0
    ln_query_size = math.log(query_size) if query_size > 0 else 0
    ln_evalue = ln_total_db_size + ln_query_size - (B * ln2)

    # Gestion des valeurs extrêmes de ln_evalue
    if ln_evalue > 700:
        evalue = float('inf')  # e-value très grande
    elif ln_evalue < -700:
        evalue = 0.0  # e-value extrêmement petite
    else:
        evalue = math.exp(ln_evalue)

    #evalue = total_db_size * query_size * (2 ** (-B))

    # Vérifier si l'HSP est significatif
    significant = evalue <= ss

    # Ajouter les informations au HSP
    hsp["bitscore"] = B
    hsp["evalue"] = evalue
    hsp["significant"] = significant

    return hsp
"""

def calculate_bitscore_and_evalue(hsp, total_db_size, query_size, ss=1e-3):
    """
    Calcule le bitscore et la e-value d'un HSP.

    Args:
        hsp (dict): Un HSP avec les champs "score" et autres.
        total_db_size (int): Taille totale des séquences dans la base de données.
        query_size (int): Taille de la séquence cible.
        ss (float): Seuil significatif pour la e-value.

    Returns:
        dict: HSP enrichi avec les champs "bitscore", "evalue", et "significant".
    """
    S = hsp["score"]  # Score brut
    K = 0.176  # Constante
    lambda_ = 0.192  # Constante
    ln2 = math.log(2)  # ln(2)

    # Calcul du bitscore
    B = (lambda_ * S - math.log(K)) / ln2
    B = max(round(B), 0)  # Bitscore ne peut pas être négatif

    # Calcul de la e-value
    evalue = total_db_size * query_size * (2 ** -B)

    # Vérifier si l'HSP est significatif
    significant = evalue <= ss

    # Ajouter les informations au HSP
    hsp["bitscore"] = B
    hsp["evalue"] = evalue
    hsp["significant"] = significant

    return hsp


#9- On trie les hsp
"""Vieille Version
def process_significant_hsps(row, total_db_size, ss=1e-3):
    query_size = len(row["Sequence"])
    fused_hsps = row["FusedHSP"]

    # Calculer bitscore et e-value pour chaque HSP
    enriched_hsps = [
        calculate_bitscore_and_evalue(hsp, total_db_size, query_size, ss)
        for hsp in fused_hsps
    ]

    # Filtrer les HSPs significatifs
    significant_hsps = [h for h in enriched_hsps if h["significant"]]

    # Ne garder que le meilleur HSP significatif par bitscore
    if significant_hsps:
        best_hsp = max(significant_hsps, key=lambda h: h["bitscore"])
        return best_hsp
    else:
        return None
"""

def process_significant_hsps(row, total_db_size, ss=1e-3):
    query_size = len(row["Sequence"])
    fused_hsps = row["FusedHSP"]

    # Calculer bitscore et e-value pour chaque HSP
    enriched_hsps = [
        calculate_bitscore_and_evalue(hsp, total_db_size, query_size, ss)
        for hsp in fused_hsps
    ]

    # Filtrer les HSPs significatifs
    significant_hsps = [h for h in enriched_hsps if h["significant"]]

    # Ne garder que le meilleur HSP significatif par bitscore
    if significant_hsps:
        best_hsp = max(significant_hsps, key=lambda h: h["bitscore"])
        return best_hsp
    else:
        # Valeur par défaut
        return {"query_start": 0, "query_end": 0, "db_start": 0, "db_end": 0,
                "aligned_query": "", "aligned_db": "", "score": 0,
                "bitscore": 0, "evalue": float("inf"), "significant": False}


#DÉFINITIONS :
# S = séquences
# g = graine
# k = taille des mots
# |g| = k
#E = pénalité

#S appartient à u_seq
g = "11111111111"
k = len(g)
E = 4

""" 
Section pour utiliser les fonctions
"""
#1-On charge les fichiers fasta
df_dataBase = read_fasta("tRNAs.fasta")
u_seq = read_fasta("unknown.fasta")

total_db_size = df_dataBase["Sequence"].apply(len).sum()


#2-On trouve les kmers
#Pour toutes les séquences, on ajoute leurs kmers à la df
#On génère des mots de taille k
"""Vieille Version
u_seq["Kmers"] = u_seq["Sequence"].apply(
    lambda seq: generate_kmers(seq, k)
    )

#3-On trouve les kmers dans la db.
u_seq["Matches"] = u_seq["Kmers"].apply(
    lambda kmers: find_kmers_in_db(kmers, df_dataBase)
    )
"""
u_seq["Kmers"] = u_seq["Sequence"].apply(lambda seq: generate_kmers_with_seed(seq, g))

# Recherche des k-mers dans la base
u_seq["Matches"] = u_seq.apply(
    lambda row: find_kmers_with_seed(row["Kmers"], df_dataBase, g), axis=1
)


#Les matchs sont enregistrés de la facon suivante :
""" 
EXEMPLE: ('GGTTCGAATCC', 52, 'M|cat|Mesostigma_viride', 52),
1- Le k-mer trouvé (GGTTCGAATCC).
2- La position du k-mer dans la séquence cible (52).
3- L'identifiant complet de la séquence dans la base (M|cat|Mesostigma_viride ou M|cat|Mesostigma_viride_2).
4- La position dans la séquence de la base (52, 53, etc.).
"""

#4/5-On extend les HSP
u_seq["Extensions"] = u_seq.apply(
    lambda row: process_hsp_extensions(row, df_dataBase, E=4), axis=1
    )


#6/7-On merge les HSP
#Mais, on doit convertir avant.
u_seq["FusedHSP"] = u_seq["Extensions"].apply(
    lambda extensions: merge_overlapping_hsps(prepare_hsps(extensions))
)

#8/9-On ajoute le best HSP au df
u_seq["BestHSP"] = u_seq.apply(
    lambda row: process_significant_hsps(row, total_db_size, ss=1e-3),
    axis=1
)


"""SECTION VISUALISTION : """
#RAPPEL : voici comment la df des S (u_seq) est construite : 
#ID, Anticodon, Info, Sequence, Kmers, Matches, Extensions, FusedHSP, BestHSP
#1-Print
#Boucle de vérification qui passe à travers chaque rangées (aka change S)
for i in range(len(u_seq)):
    print("Pour la séquence S : " + u_seq.loc[i, 'Sequence'])
    print("On a les caractéristiques suivantes : ")
    print("Kmers : \n" + str(u_seq.loc[i, 'Kmers']))
    print("Matches : \n" + str(u_seq.loc[i, 'Matches']))
    print("Extensions : \n" + str(u_seq.loc[i, 'Extensions']))
    print("FusedHSP : \n" + str(u_seq.loc[i, 'FusedHSP']))
    print("BestHSP : \n" + str(u_seq.loc[i, 'BestHSP']))
    
#2-Fichiers
#Chaque fois qu'on modifie la df, on génère un csv pour la visualiser
u_seq.to_csv("repr_DF", index=False)
df_dataBase.to_csv("df_DataBase", index=False)

"""Vérification rapide de merge_overlapping_hsps"""
# hsps = [
#     {"query_start": 10, "query_end": 20, "db_start": 15, "db_end": 25, 
#      "aligned_query": "ACGTACGTAC", "aligned_db": "ACGTACGTAC", "score": 50},
#     {"query_start": 18, "query_end": 30, "db_start": 23, "db_end": 35, 
#      "aligned_query": "TACGTTACGTA", "aligned_db": "TACGTTACGTA", "score": 60},
# ]
# merged = merge_overlapping_hsps(hsps)
# print(merged)
