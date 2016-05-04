wget 
unzip Example.zip
export PYTHONPATH=`pwd`
./qc.py -s OBP



0 - Installation
----------------
Pré-requis :
- python (version 2.xx; pas testé sous python 3)
- obspy 1.0.1 (www.obspy.org)

Déclarer la variable PYTHONPATH : 
export PYTHONPATH=`pwd`
Ou la mettre dans votre .bashrc

1- Créer une nouvelle station dans le fichier station_dictionnary.py
--------------------------------------------------------------------
C'est un dictionnaire python qui regroupe différentes informations utiles aux codes pour chaque station. 
Le mieux est de faire un copier/coller d'une station.
Les clés obligatoires (c'est à dire où l'information doit être correcte) sont :
- network : code réseau présent dans les fichiers mseed
- locid : location code présent dans les fichiers mseed
- channels : nom des channels devant être traités (et présents dans les fichiers mseed)
- path_data : chemin générique vers les fichiers de données mseed (cf 2.2- chemin vers les données)

2- Préparation des données
--------------------------
2.1 - Fichiers journaliers

Le code nécessite des données au format mseed (ça doit fonctionner avec d'autres formats) avec un fichier par jour démarrant, si possible, à une heure entière (c'est pas hyper critique).
Si vous avez déjà ça : passer à l'étape 2.2
Sinon, vous pouvez utiliser le code "extract2sds.py" pour fabriquer des fichiers journaliers à partir de données stockées dans un répertoire et stocker le tout dans une structure de fichier "SDS". C'est surtout utile pour rapidement ré-organiser des données "brutes" produites par un Q330 (après dé-archivage via sdrsplit par exemple).
C'est un clone (très simplifié!) de qmerge. Il évidemment est préférable d'utiliser qmere ou msmod.
Exemple :
python extract2sds.py -f "/home/toto/RAW/NEUF/QT*" -s NEUF -n XX -l 00 -o "/home/toto/Data/"
va :
-lire tous les fichiers démarrant par QT dans le répertoire /home/toto/RAW/NEUF
-extraire des fichiers journaliers de données et les mettre dans /home/toto/Data/201?/XX/NEUF/HH?.D/
-renommer (dans les header miniseed) la station en NEUF, le réseau en XX et le location code en 00
Possibilité de restreindre l'extraction à une période de temps avec les options -b et -e (example -b 20120312 -e 20120315 extrait uniquement entre le 12 et 15 mars 2012)

2.2 Déclaration du chemin vers les données mseed

Dans le fichier station_dictionnary.py, remplir la clé 'path_data' avec le bon chemin et en utilisant les variables %(sta), %(chan), ... pour remplacer automatiquement les noms de station, channel, ...
Exemple : Mes données sont dans un répertoire du type /home/toto/Data/2012/XX/NEUF/HHZ.D/NEUF.XX.HHZ.00.D.2012.214
Mettre dans 'path_data' : /home/toto/Data/%(year)s/%(net)s/%(sta)s/%(chan)s.D/*%(year)s.%(day)s.
!!! attention au "s" (obligatoire)

2.3 Réponse instrumentale

2.3.1 J'ai un fichier dataless
C'est l'idéal. Par contre la réponse doit y être définie en radian/s (flag "B") et non pas en Hertz (flag "A") car obspy ne sait pas (encore) lire ce type de dataless.
Indiquer simplement le chemin complet du dataless pour la clé 'dataless_file' de station_dictionnary.py 

2.3.2 J'ai pas de fichier dataless
C'est mal ... mais pas critique ! car vous connaissez au moins le type de numériseur et de capteur utilisé (si pas esotérique).
Dans ce cas remplissez les clés 'digitizer' et 'sensor' de station_dictionnary.py. Les réponses (poles/zeros/gains) correspondantes seront lues dans le dictionnaire/fichier instruments.py. Il faut donc utiliser des noms de 'digitizer' et 'sensor' qui existent dans instruments.py
Rq : Seul le lsb des numériseurs sera utilisé pour la déconvolution => réponse fausse proche de la fréquence de Nyquist

(A VENIR : make_dataless_from_generic.py)


3. Lancer le code
-----------------
Le code principale est run_qc_arg
Dans le terminal taper run_qc_arg -h pour avoir toutes les options
En gros ce code va :
- lire les données à partir des infos trouvées dans station_dictionnary.py
- Faire l'analyse de bruit via le module qc.py (clone de ppsd de obspy) pour des périodes de 7 jours (cf option -nb_days_pack)
- Stocker les résultats dans un fichier au format .pkl (pickle) dans le répertoire PATH_PKL
- Fabriquer la figure de synthèse, au format png, dans le répertoire PATH_PLT
Les noms des fichiers .pkl et .plt sont de la forme [network].[station].[location].[channel].{pkl,png}
Les répertoires par défault pour ces fichiers (PATH_PKL et PATH_PLT) sont définis dans le fichier default_qc_path.py. Ils peuvent être forcés via les options -pkl et -plt

4. Re-génération des figures
----------------------------
si vous voulez refaire les figures pour une période de temps donnée, avec ou sans la réprésentation des classes de sites, ... vous pouvez utiliser plot_qc.py
Cela évite de retraiter (relire) toutes les données. Ce code va simplement lire un ensemble de fichier .pkl et utiliser les différents options de la méthode plot de laclasse QC
Exemple :
python plot_qc.py -s ST1 ST2 -c HHZ -b 20130101 -e 20130115 -pkl /home/toto/PKL -nc
Va regénérer les figures pour les composantes HHZ des stations ST1 et ST2, pendant la période allant du 01/01/2013 au 15/01/2013, sans représenter les traits des classes A et B, à partir des fichiers pkl correspondants présents dans /home/toto/PKL. Les nouvelles figures png seront stockées dans le répertoire courant (pas de vérification si le fichier existe déjà)


Example :
Suite de commandes pour traiter la station OBP (donnée dans Example) :

export PYTHONPATH=`pwd`


./extract2sds.py -f "Example/Raw_Data/ROVI/*" -s ROVI -n XX -l 00 -o "Example/MSEED/"

./run_qc_arg -s ROVI

./plot_qc.py -s ROVI -b 20120919190000 -e 2012092020000 -plt "Example/PLT/Temp/" -nc


5. Autres codes (experimental)
------------------------------
qc_comp_2sta.py : Pour comparer 2 stations (ou 2 channels) sur un laps de temps donné
qc_comp_2periods.py : Pour comparer une station/channel sur deux laps de temps donnés
qc_comp_to_ref.py : Pour effectuer un classement de stations sur un laps de temps donné en terme de niveau de bruit par rapport au NLNM dans plusieurs gammes de périodes définies. Permet également d'afficher les sites repondant aux critères des classes A et B dans ces gammes de période.

RQ : Pour ces codes, il est nécessaire de modifier la partie "PARAMETERS" dans les programmes correspondants. Ils sont exécutable dans l'environnement ipython
TODO : Passer des arguments, généraliser, vérifier, ...



modifier:
- les path dans default_qc_path.py
- les dates min et max de run_qc_arg ET extract2sds.py
- le lower de run_qc_arg

acq=eval(str(eval(sta)['digitizer']).lower)
sismo=eval(str(eval(sta)['sensor']).lower)

