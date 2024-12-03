# File Processor

## Utilisation
1. Téléchargez l'exécutable pour votre système d'exploitation.
2. Lancez-le (double-cliquez sur `app.exe` sur Windows, ou `./app` sur Linux/macOS).
3. Une page web s'ouvrira automatiquement dans votre navigateur.
4. Chargez vos fichiers `.txt` ou `.fasta` et cliquez sur "Process Files" pour afficher les résultats.

## Dépendances
Aucune, l'exécutable est autonome.  
Dans tous les cas, il y a un fichier requirements contenant les installations.  


## À NOTER
- Le projet est compilé avec pyinstaller.
Pour chaque OS : on compile avec : pyinstaller --onefile --add-data "templates:index" --add-data "static:index" app.py

## Structure
projet
|
 - static
 |
