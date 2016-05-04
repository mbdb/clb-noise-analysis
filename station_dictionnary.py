# Station dictionnary for run_qc.py and make_dataless_from_generic.py
# channels are the channels that will be processed in run_qc.py
# Each station consists in a seismometer (sensor) and a digitizer
# Sensor's and digitizer's PAZ have to be declared in the file instruments.py
# if no sensor+digitizer is given run_qc expects to find a dataless file (dataless_file full path)
# path_data define the directory where are stored the mseed files. Generic formatters are :
#	station : %(sta)s 
#	network : %(net)s 
#	locid : %(locid)s 
#	channel : %(chan)s 
#	day : %(day)s   !!! MANDATORY !!!
#	year : %(year)s


OBP = {'network' : 'G',
	'channels' : ['HHZ', 'HHN','HHE' ],  
	'locid' : '10' ,
	'digitizer' : 'Q330HR',
	'sensor' : 'STS2' ,
    'dataless_file' : 'Example/Dataless/OBP.dataless' ,
	'path_data' : 'Example/SDS/%(year)s/%(net)s/%(sta)s/%(chan)s.D/*.%(day)s' }
