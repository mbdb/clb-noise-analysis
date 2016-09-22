
1- Crer une nouvelle station dans le fichier station_dictionnary.py
--------------------------------------------------------------------
C'est un dictionnaire python qui regroupe diffrentes informations utiles aux codes pour chaque station. 
Le mieux est de faire un copier/coller d'une station.
Les cls obligatoires (c'est  dire o l'information doit tre correcte) sont :
- network : code rseau prsent dans les fichiers mseed
- locid : location code prsent dans les fichiers mseed
- channels : nom des channels devant tre traits (et prsents dans les fichiers mseed)
- path_data : chemin gnrique vers les fichiers de donnes mseed (cf 2.2- chemin vers les donnes)

2- Preparation des donnes
--------------------------
2.1 - Fichiers journaliers

Le code ncessite des donnes au format mseed (a doit fonctionner avec d'autres formats) avec un fichier par jour dmarrant, si possible,  une heure entire (c'est pas hyper critique).
Si vous avez dj a : passer  l'tape 2.2
Sinon, vous pouvez utiliser le code "extract2sds.py" pour fabriquer des fichiers journaliers  partir de donnes stockes dans un rpertoire et stocker le tout dans une structure de fichier "SDS". C'est surtout utile pour rapidement r-organiser des donnes "brutes" produites par un Q330 (aprs d-archivage via sdrsplit par exemple).
C'est un clone (trs simplifi!) de qmerge. Il videmment est prfrable d'utiliser qmere ou msmod.
Exemple :
python extract2sds.py -f "/home/toto/RAW/NEUF/QT*" -s NEUF -n XX -l 00 -o "/home/toto/Data/"
va :
-lire tous les fichiers dmarrant par QT dans le rpertoire /home/toto/RAW/NEUF
-extraire des fichiers journaliers de donnes et les mettre dans /home/toto/Data/201?/XX/NEUF/HH?.D/
-renommer (dans les header miniseed) la station en NEUF, le rseau en XX et le location code en 00
Possibilite de restreindre l'extraction  une priode de temps avec les options -b et -e (example -b 20120312 -e 20120315 extrait uniquement entre le 12 et 15 mars 2012)

2.2 Declaration du chemin vers les donnes mseed

Dans le fichier station_dictionnary.py, remplir la cle 'path_data' avec le bon chemin 
2.3 Reponse instrumentale

2.3.1 J'ai un fichier dataless
C'est l'idal. Par contre la reponse doit y tre definie en radian/s (flag "B") et non pas en Hertz (flag "A") car obspy ne sait pas (encore) lire ce type de dataless.
Indiquer simplement le chemin complet du dataless pour la cle 'dataless_file' de station_dictionnary.py 

2.3.2 J'ai pas de fichier dataless
C'est mal ... mais pas critique ! car vous connaissez au moins le type de numeriseur et de capteur utilise (si pas esotrique).
Dans ce cas remplissez les cles 'digitizer' et 'sensor' de station_dictionnary.py. Les reponses (poles/zeros/gains) correspondantes seront lues dans le dictionnaire/fichier instruments.py.
Il faut donc utiliser des noms de 'digitizer' et 'sensor' qui existent dans instruments.py
Rq : Seul le lsb des numeriseurs sera utilise pour la deconvolution => reponse fausse proche de la frequence de Nyquist

(A VENIR : make_dataless_from_generic.py)


3. Lancer le code
-----------------
Le code principale est run_qc_arg
Dans le terminal taper run_qc_arg -h pour avoir toutes les options
En gros ce code va :
- lire les donnes  partir des infos trouves dans station_dictionnary.py
- Faire l'analyse de bruit via le module qc.py (clone de ppsd de obspy) pour des priodes de 7 jours (cf option -nb_days_pack)
- Stocker les rsultats dans un fichier au format .pkl (pickle) dans le rpertoire PATH_PKL
- Fabriquer la figure de synthse, au format png, dans le rpertoire PATH_PLT
Les noms des fichiers .pkl et .plt sont de la forme [network].[station].[location].[channel].{pkl,png}
Les rpertoires par dfault pour ces fichiers (PATH_PKL et PATH_PLT) sont dfinis dans le fichier default_qc_path.py. Ils peuvent tre forcs via les options -pkl et -plt

4. Re-genration des figures
----------------------------
si vous voulez refaire les figures pour une priode de temps donne, avec ou sans la rprsentation des classes de sites, ... vous pouvez utiliser plot_qc.py
Cela vite de retraiter (relire) toutes les donnes. Ce code va simplement lire un ensemble de fichier .pkl et utiliser les diffrents options de la mthode plot de laclasse QC
Exemple :
python plot_qc.py -s ST1 ST2 -c HHZ -b 20130101 -e 20130115 -pkl /home/toto/PKL -nc
Va regnrer les figures pour les composantes HHZ des stations ST1 et ST2, pendant la priode allant du 01/01/2013 au 15/01/2013, sans reprsenter les traits des classes A et B,  partir des fichiers pkl correspondants prsents dans /home/toto/PKL. Les nouvelles figures png seront stockes dans le rpertoire courant (pas de vrification si le fichier existe dj)


Example :
Suite de commandes pour traiter la station OBP (donne dans Example) :

export PYTHONPATH=`pwd`


./extract2sds.py -f "Example/Raw_Data/ROVI/*" -s ROVI -n XX -l 00 -o "Example/MSEED/"

./run_qc_arg -s ROVI

./plot_qc.py -s ROVI -b 20120919190000 -e 2012092020000 -plt "Example/PLT/Temp/" -nc


5. Autres codes (experimental)
------------------------------
qc_comp_2sta.py : Pour comparer 2 stations (ou 2 channels) sur un laps de temps donn
qc_comp_2periods.py : Pour comparer une station/channel sur deux laps de temps donns
qc_comp_to_ref.py : Pour effectuer un classement de stations sur un laps de temps donn en terme de niveau de bruit par rapport au NLNM dans plusieurs gammes de priodes dfinies. Permet galement d'afficher les sites repondant aux critres des classes A et B dans ces gammes de priode.

RQ : Pour ces codes, il est ncessaire de modifier la partie "PARAMETERS" dans les programmes correspondants. Ils sont excutable dans l'environnement ipython
TODO : Passer des arguments, gnraliser, vrifier, ...



modifier:
- les path dans default_qc_path.py
- les dates min et max de run_qc_arg ET extract2sds.py
- le lower de run_qc_arg

acq=eval(str(eval(sta)['digitizer']).lower)
sismo=eval(str(eval(sta)['sensor']).lower)


./qc.py  -s BUFF10    -pkl ~/git/clb-noise-analysis/PKL/ -plt ~/git/clb-noise-analysis/PLT/  -b 2016201

./qc.py  -s MEUD00 BUFF10   


./plot_qc.py  -s BUFF -c BHZ -n XX -l 00 -pkl ~/git/clb-noise-analysis/PKL/ -plt ~/git/clb-noise-analysis/PLT/  -b 2016-07-23
