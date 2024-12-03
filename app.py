import logging
from flask import Flask, render_template, request, jsonify
import pandas as pd
from TP3_3295 import *

# Configuration des journaux
logging.basicConfig(level=logging.INFO)

# Initialisation de l'application Flask
app = Flask(__name__)

#Fonction principale pour traiter les fichiers
def process_files(query_path, db_path, g, E, ss):
    logging.info("Lecture des fichiers FASTA")
    df_dataBase = read_fasta(db_path)
    u_seq = read_fasta(query_path)

    # Calcul des paramètres
    logging.info("Calcul des paramètres")
    total_db_size = df_dataBase["Sequence"].apply(len).sum()
    k = len(g)

    # Étape 1 : Générer les kmers
    logging.info("Génération des kmers")
    u_seq["Kmers"] = u_seq["Sequence"].apply(lambda seq: generate_kmers_with_seed(seq, k))

    # Étape 2 : Trouver les kmers dans la base
    logging.info("Recherche des kmers dans la base")
    u_seq["Matches"] = u_seq["Kmers"].apply(lambda kmers: find_kmers_with_seed(kmers, df_dataBase))

    # Étape 3 : Étendre les HSPs
    logging.info("Extension des HSPs")
    u_seq["Extensions"] = u_seq.apply(lambda row: process_hsp_extensions(row, df_dataBase, E=E), axis=1)

    # Étape 4 : Fusionner les HSPs
    logging.info("Fusion des HSPs")
    u_seq["FusedHSP"] = u_seq["Extensions"].apply(lambda extensions: merge_overlapping_hsps(prepare_hsps(extensions)))

    # Étape 5 : Identifier le meilleur HSP significatif
    logging.info("Identification des meilleurs HSPs significatifs")
    u_seq["BestHSP"] = u_seq.apply(lambda row: process_significant_hsps(row, total_db_size, ss=ss), axis=1)

    logging.info("Traitement terminé")
    return u_seq.to_dict(orient="records")


# Route principale
@app.route("/")
def index():
    return render_template("index.html")

# Route pour traiter les fichiers
@app.route("/process", methods=["POST"])
def process():
    try:
        logging.info("Début du traitement")

        # Charger les fichiers envoyés
        query_file = request.files["query_file"]
        db_file = request.files["db_file"]

        # Sauvegarder temporairement les fichiers pour traitement
        query_path = "temp_query.fasta"
        db_path = "temp_db.fasta"
        query_file.save(query_path)
        db_file.save(db_path)

        # Récupérer les nouvelles variables
        g = request.form.get("g", "11111111111")  # Graine par défaut
        E = int(request.form.get("E", 4))        # Pénalité par défaut
        ss = float(request.form.get("ss", 0.001))  # Seuil par défaut

        # Appeler la fonction principale
        u_seq = process_files(query_path, db_path, g, E, ss)

        # Convertir la DataFrame en HTML
        u_seq_html = pd.DataFrame(u_seq).to_html(index=False)

        # Retourner le HTML au frontend
        return jsonify({"table": u_seq_html})
    except Exception as e:
        logging.error(f"Erreur lors du traitement : {e}")
        return jsonify({"error": str(e)}), 500




if __name__ == "__main__":
    app.run(debug=True)
