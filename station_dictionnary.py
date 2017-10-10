# Station dictionnary for run_qc.py and make_dataless_from_generic.py
# each station name should be unique like OBP10 and OBP00
# channels are the channels that will be processed in run_qc.py
# Each station consists in a seismometer (sensor) and a digitizer or a dataless_fill
# Sensor's and digitizer's PAZ have to be declared in the file instruments.py
# if no sensor+digitizer is given run_qc expects to find a dataless file (dataless_file full path)
# path_data define the directory where are stored the mseed files in a SDS archive

OBP10 = {'network': 'G',
         'station':'OBP',
        'locid': '10',
        'channels': ['HHZ', 'HHN', 'HHE'],
        'digitizer': 'q330hr',
        'sensor': 'sts2',
        'path_data': 'Example/SDS'}


MEUD00 = { 'network': 'XX',
          'station':'MEUD',
          'locid': '00',
          'channels': ['HLZ', 'HLN', 'HLE'],
          'dataless_file': 'Example/Dataless/MEUD_XX.dataless',
          'path_data': 'Example/SDS',
          'title_comment': r'cloche $\mu$-metal'
}
