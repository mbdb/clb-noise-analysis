# clb-noise-analysis

Requirements:
+ OBSPY 1.0.1

1) Each station to be analysed needs to be in the station_dictionnary.py file with a unique name. Just copy-paste the example and modify
- station
- network
- locid
- channels to be analysed
- dataless_file or digitizer and sensor (instruments.py)
- path_data of the SDS structure where data is
- optional title_comment

2) Daily data files needs to be in a SDS structure. If you don't have it already, you can use extract2sds.py to create it. This programm can at the same time change the station name, the locid, the network in the files. 
Ex:
``` code
./extract2sds.py -f "Example/Raw_Data/VINC/*" -s VINC -n XX -l 00 -o Example/SDS/
````
For each station you need a dataless file which could be created with PDCC (https://ds.iris.edu/ds/nodes/dmc/software/downloads/pdcc)
or a sensor and a digitizer in the instruments.py file

3) To analyse the data:
``` code
./qc.py -s MEUD00
````
You can analyse multiple stations
``` code
./qc.py -s MEUD00 OBP10
````

You can analyse with specific dates
``` code
./qc.py -s OBP10 -b 2015,261 -e 2016,264
````


Example:
You can download some files to test:
``` code
wget http://www.ipgp.fr/~bonaime/Example.zip
unzip Example.zip



```

![image](doc/XX.MEUD.00.HLZ.png)

4) To plot a specific period of the data, you can use plot_qc.py class:
``` code
 ./plot_qc.py  -s BUFF -c BHZ -n XX -l 00 -b 2016-05-29 -e 2016-06-29
``` 
