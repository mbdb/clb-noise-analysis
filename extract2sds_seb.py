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
from obspy.core import *
from datetime import datetime
from default_qc_path import *
import tempfile
import sys

# ----------- ARGUMENTS  / PARAMETERS --------------------------
CHAN_default = ['HHZ', 'HHN', 'HHE', 'BHZ', 'BHN', 'BHE', 'LHZ', 'LHN', 'LHE']

argu_parser = argparse.ArgumentParser(
    description='extract one hour long segments from FILES and store them in the PATH_SDS directory. Use options -s -n and -l to force the station, network and location code in mseed files. Otherwise the orginal informations are preserved. Directory''s name of the SDS is based on these codes. ')
argu_parser.add_argument("-f", "--files", required=True,
                         help="File(s) to be read. If it contains wildcards(*,?,...), the string must be quoted. Example : python extract2sds.py -s TOTO -n XX \"../*HH?.*\"")
argu_parser.add_argument("-c", "--channels", nargs='+', default=CHAN_default,
                         help="Process only CHANNELS. Do not use this option if you want to process all available channels for the corresponding station. Separate CHANNELS with spaces. No wildcard. Default is all LH?,BH?,HH?")
argu_parser.add_argument("-b", "--starttime", default=UTCDateTime(2000, 1, 1), type=UTCDateTime,
                         help="Start time for processing. Various format accepted. Example : 2012,2,1 / 2012-02-01 / 2012,032 / 2012032 / etc ... See UTCDateTime for a complete list. Default is 2010-1-1")
argu_parser.add_argument("-e", "--endtime", default=UTCDateTime(2099, 1, 1), type=UTCDateTime,
                         help="End time for processing. Various format accepted. Example : 2012,2,1 / 2012-02-01 / 2012,032 / 2012032 / etc ... See UTCDateTime for a complete list. Default is 2015-1-1")
argu_parser.add_argument("-s", "--station", help="set station code")
argu_parser.add_argument("-n", "--network", help="set network code")
argu_parser.add_argument(
    "-l", "--location", help="set location code. \"\" if none ")
argu_parser.add_argument("-o", "--path_sds", default=PATH_SDS, help="Base directory for the SDS. Default is " +
                         PATH_SDS + " (default path can be modified in default_qc_path)")

args = argu_parser.parse_args()
start = args.starttime
stop = args.endtime
PATH_SDS = args.path_sds + "/"
files = args.files
CHAN = args.channels
sta = args.station
net = args.network
locid = args.location

# ----------- END OF ARGUMENTS / PARAMETERS -------------------


# Create tmp directory to store links to valid mseed files
path_link = tempfile.mkdtemp()
file_list = []

# Read header information for selected streams

for f in glob.glob(files):
    try:
        Temp = read(f, headonly=True)
    except:
        #       pass
        print "unable to read " + f
    else:
        chan = np.unique([s.stats['channel'] for s in Temp])
        if len(np.intersect1d(chan, CHAN)) > 0:
            file_list.append(f)
        else:
            print "No valid CHANNEL in " + os.path.basename(f)
Files = path_link + '/*'
try:
    Stream_in_head = Stream()
    for file in file_list:
        Stream_in_head += read(file, headonly=True)
except:
    e = sys.exc_info()[0]
    print e
    print("No valid data file found for " + files)
else:
    print Stream_in_head
    chan_exist = np.unique([s.stats['channel'] for s in Stream_in_head])
    Start_tr = [s.stats['starttime'] for s in Stream_in_head]
    Stop_tr = [s.stats['endtime'] for s in Stream_in_head]
    net_dir = np.unique([s.stats['network'] for s in Stream_in_head])
    sta_dir = np.unique([s.stats['station'] for s in Stream_in_head])
    CHAN = np.intersect1d(chan_exist, CHAN)
    tmin = min(Start_tr)
    tmax = max(Stop_tr)
    year_dir = np.arange(tmin.year, tmax.year + 1)
    if net:
        net_dir = np.array([net])
    if sta:
        sta_dir = np.array([sta])
    # Create output directories (SDS structure)
    for y in year_dir:
        try:
            os.mkdir(PATH_SDS + str(y))
        except:
            pass
        for n in net_dir:
            try:
                os.mkdir(PATH_SDS + str(y) + '/' + n)
            except:
                pass
            for s in sta_dir:
                try:
                    os.mkdir(PATH_SDS + str(y) + '/' + n + '/' + s)
                except:
                    pass
                for chan in CHAN:
                    try:
                        os.mkdir(PATH_SDS + str(y) + '/' + n +
                                 '/' + s + '/' + chan + ".D")
                    except:
                        pass

    # Create array of Days encompassed
    tmin = UTCDateTime(tmin.timestamp // 86400 * 86400)
    tmax = UTCDateTime(tmax.timestamp // 86400 * 86400 + 86400)
    Days = np.arange(max(start, tmin), min(tmax, stop), 86400)

    # Process every days for each channel
    for day in Days:
        Stream = read(Files, starttime=day, endtime=day +
                      86400, nearest_sample=False)
        for chan in CHAN:
            S = Stream.select(channel=chan)
            # Merge
            # S.merge()
            S.merge(fill_value="latest")
            # change ID info
            for tr in S:
                # fill masked array
                # if np.ma.isMaskedArray(tr.data):
                #   tr.data=tr.data.filled()
                if sta:
                    tr.stats.station = sta
                else:
                    sta = tr.stats.station
                if net:
                    tr.stats.network = net
                else:
                    net = tr.stats.network
                if locid:
                    tr.stats.location = locid
                else:
                    locid = tr.stats.location
                t0 = tr.stats.starttime
                nameout = "%s/%04i/%s/%s/%s.D/%s.%s.%s.%s.D.%04i.%03i" % (
                    PATH_SDS, t0.year, net, sta, chan, sta, net, chan, locid, t0.year, t0.julday)
                tr.write(nameout, format="MSEED")
                print("write %s (%i samples, %10.2f seconds)" %
                      (nameout, tr.stats.npts, tr.stats.endtime - t0))
