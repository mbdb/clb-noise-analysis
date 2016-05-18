#!/usr/bin/env python
#------------------------------------------------------------------------------
# Filename: extract2sds.py
#  Version: 1.2
#  Purpose: Convert seismic data into 1 day long segments starting at time 00h00min00s and store them into a Standard Data Structure (SDS)
#     Note: Input data must be in format readable with obspy (works for mseed,sac, ... as well as raw data from Quanterra digitizers)
#           Unreadable files are ommitted
#           Temporary files (links) are written to /tmp. Thus this directory must exist
#    Author: Jerome Vergne
#    Email: jerome.vergne@unsitra.fr
#
#------------------------------------------------------------------------------
import os
import glob
import argparse
import numpy as np
from obspy.core import read
from obspy import Stream
from obspy import UTCDateTime
from default_qc_path import PATH_SDS


# ----------- ARGUMENTS  / PARAMETERS --------------------------
CHAN_default = ['HHZ', 'HHN', 'HHE', 'BHZ', 'BHN', 'BHE', 'LHZ', 'LHN', 'LHE']

argu_parser = argparse.ArgumentParser(
    description='extract one hour long segments from FILES and store them in\
    the PATH_SDS directory. Use options -s -n and -l to force the station, \
    network and location code in mseed files. Otherwise the orginal \
    informations are preserved. Directory''s name of the SDS is based \
    on these codes. ')
argu_parser.add_argument("-f", "--files", required=True,
                         help="File(s) to be read. If it contains wildcards\
                         (*,?,...), the string must be quoted. Example : python\
                          extract2sds.py -s TOTO -n XX \"../*HH?.*\"")
argu_parser.add_argument("-c", "--channels", nargs='+', default=CHAN_default,
                         help="Process only CHANNELS. Do not use this option if\
                          you want to process all available channels for the\
                           corresponding station. Separate CHANNELS with spaces\
                           . No wildcard. Default is all LH?,BH?,HH?")
argu_parser.add_argument("-b", "--starttime", default=UTCDateTime(2000, 1, 1), type=UTCDateTime,
                         help="Start time for processing. Various format \
                         accepted. Example : 2012,2,1 / 2012-02-01 / 2012,032\
                          / 2012032 / etc ... See UTCDateTime for a complete\
                           list. Default is 2010-1-1")
argu_parser.add_argument("-e", "--endtime", default=UTCDateTime(2099, 1, 1), type=UTCDateTime,
                         help="End time for processing. Various format accepted\
                         . Example : 2012,2,1 / 2012-02-01 / 2012,032 / 2012032\
                          / etc ... See UTCDateTime for a complete list.\
                           Default is 2015-1-1")
argu_parser.add_argument("-s", "--station", help="set station code")
argu_parser.add_argument("-n", "--network", help="set network code")
argu_parser.add_argument(
    "-l", "--location", help="set location code. \"\" if none ")
argu_parser.add_argument("-o", "--path_sds", default=PATH_SDS, help="Base\
 directory for the SDS. Default is " +PATH_SDS + " (default path can be \
 modified in default_qc_path)")

args = argu_parser.parse_args()
start = args.starttime
stop = args.endtime
PATH_SDS = args.path_sds + "/"
input_files = args.files
CHAN = args.channels

sta = args.station
net = args.network
locid = args.location

# ----------- END OF ARGUMENTS / PARAMETERS -------------------

# read all files
all_streams = Stream()
for input_file in glob.glob(input_files):
    print "Reading data from " + input_file
    all_streams += read(input_file, nearest_sample=False,
                        sourcename='*.[BHL][HL][ZNE]', details=True)

# Create array of Days encompassed
min_startime = UTCDateTime(
    (np.min([f.stats.starttime for f in all_streams.traces])).timestamp // 86400 * 86400)
max_startime = UTCDateTime((np.max(
    [f.stats.starttime for f in all_streams.traces])).timestamp // 86400 * 86400 + 86400)
Days = np.arange(max(start, min_startime), min(max_startime, stop), 86400)

# Process every days for each channel
for day in Days:
    day_streams = all_streams.slice(day, day + 86400)

    for day_channel in CHAN:
        day_stream = day_streams.select(channel=day_channel)
        day_stream.merge(fill_value="latest")
        # change ID info
        for day_trace in day_stream:
            if sta:
                day_trace.stats.station = sta
            else:
                sta = day_trace.stats.station
            if net:
                day_trace.stats.network = net
            else:
                net = day_trace.stats.network
            if locid:
                day_trace.stats.location = locid
            else:
                locid = day_trace.stats.location
            t0 = day_trace.stats.starttime
            file_dir = PATH_SDS + str(t0.year) + '/' + \
                net + '/' + sta + '/' + day_channel + '.D'
            if not os.path.exists(file_dir):
                os.makedirs(file_dir)
            nameout = "%s/%s.%s.%s.%s.D.%04i.%03i" % (
                file_dir, sta, net, day_channel, locid, t0.year, t0.julday)

            print nameout
            print day_trace
            day_trace.write(nameout, format="MSEED")
            print("write %s (%i samples, %10.2f seconds)" %
                  (nameout, day_trace.stats.npts, day_trace.stats.endtime - t0))