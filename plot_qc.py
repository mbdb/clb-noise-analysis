#!/usr/bin/env python
#------------------------------------------------------------------------------
# Filename: plot_qc.py
#  Version: 1.0
#  Purpose: small program to (re)plot the quality sheet based on pkl files computed from run_qc.py. Save the figure to a png file. This program only pass various options of the QC.plot method at command line
#  CAUTION : .pkl files must be named in the form [network].[station].[location].[channel].pkl
#   Author: Jerome Vergne
#    Email: jerome.vergne@unsitra.fr
#
#------------------------------------------------------------------------------
import os
import glob
import argparse
from qc import *
from obspy.imaging.cm import pqlx, get_cmap

import sys
PYTHON_VERSION = [int(n) for n in sys.version[:3].split('.')]
if PYTHON_VERSION < [3, 0, 0]:
    import cPickle as Pickle
else:
    import _pickle as Pickle


# PARAMETERS
# ----------

# Arguments
argu_parser = argparse.ArgumentParser(
    description="Plot quality sheets for a set of STATIONS. A corresponding .pkl file must exsit in the PKL directory")
argu_parser.add_argument("-s", "--stations", nargs='+', required=True,
                         help="list of stations name. Separate STATIONS with spaces")
argu_parser.add_argument("-c", "--channels", nargs='+',
                         help="Process only CHANNELS. Do not use this option if you want to process all available channels for the corresponding station. Separate CHANNELS with spaces. No wildcard. Default is all channels")
argu_parser.add_argument("-n", "--networks", nargs='+',
                         help="Process only NETWORKS. Do not use this option if you want to process all available networks for the corresponding station. Separate NETWORKS with spaces. No wildcard. Default is all networks")
argu_parser.add_argument("-l", "--locations", nargs='+',
                         help="Process only LOCATIONS id. Do not use this option if you want to process all available locations id for the corresponding station. Separate LOCATIONS with spaces. No wildcard. Default is all locations")
argu_parser.add_argument("-b", "--starttime", default=UTCDateTime(2001, 1, 1), type=UTCDateTime,
                         help="Start time for processing. Various format accepted. Example : 2012,2,1 / 2012-02-01 / 2012,032 / 2012032 / etc ... See UTCDateTime for a complete list. Default is 2001-1-1")
argu_parser.add_argument("-e", "--endtime", default=UTCDateTime(2050, 1, 1), type=UTCDateTime,
                         help="End time for processing. Various format accepted. Example : 2012,2,1 / 2012-02-01 / 2012,032 / 2012032 / etc ... See UTCDateTime for a complete list. Default is 2015-1-1")
argu_parser.add_argument("-pkl", "--path_pkl", default='./PKL',
                         help="Directory for pkl files. Default is ./PKL (defined in default_path_qc.py file)")
argu_parser.add_argument("-plt", "--path_plt", default="./PLT",
                         help="Output directory for plt files. Default is ./PLT (directory where you stand). CAREFULL : will overwrite")
argu_parser.add_argument("-np", "--no_show_percentiles", default=False,
                         action='store_true', help="Do not plot percentiles. Default = plot")
argu_parser.add_argument("-nc", "--no_show_class_models", default=False,
                         action='store_true', help="Do not plot A and B class limits. Default = plot")
argu_parser.add_argument("--color_map", default='pqlx',
                         help="Color map for PPSD. Default is pqlx")



args = argu_parser.parse_args()

# List of stations
stations = args.stations
channels = args.channels
locations = args.locations
networks = args.networks
# List of networks
# List of location code
# Time span
start = args.starttime
stop = args.endtime
# INPUT Path
PATH_PKL = args.path_pkl + "/"
# Output Path
PATH_PLT = args.path_plt + "/"
# Plot parameters
np = args.no_show_percentiles
nc = args.no_show_class_models


# Color Map
if args.color_map == 'pqlx':
    cmap = pqlx
elif args.color_map == 'viridis_white':
    from obspy.imaging.cm import viridis_white
    cmap = viridis_white
elif args.color_map == 'viridis_white_r':
    from obspy.imaging.cm import viridis_white_r
    cmap = viridis_white_r
else:
    cmap = get_cmap(args.color_map)

# --------- MAIN ------------
PKL_FILES = glob.glob(PATH_PKL + "*.pkl")
# Loop over files
for pkl_file in PKL_FILES:
    try:
        (net, sta, loc, chan, temp) = pkl_file.split("/")[-1].split(".")
    except:
        pass
    else:
        if isinstance(stations, list):
            test_sta = (sta in stations)
        else:
            test_sta = True
        if isinstance(channels, list):
            test_chan = (chan in channels)
        else:
            test_chan = True
        if isinstance(networks, list):
            test_net = (net in networks)
        else:
            test_net = True
        if isinstance(locations, list):
            test_loc = (loc in locations)
        else:
            test_loc = True
        if test_sta and test_chan and test_net and test_loc:
            try:
                S = Pickle.load(open(pkl_file, 'rb'))
            except:
                print("Can't read " + pkl_file)
            else:
                print("Load the pickle file : " + pkl_file)
                tb = max(start, min(S.times_used))
                te = min(stop, max(S.times_used))

                filename_plt = PATH_PLT + net + "." + sta + "." + loc + "." + chan + ".png"
                S.plot(cmap, filename=filename_plt, show_percentiles=not(np),
                       show_class_models=not(nc), starttime=tb, endtime=te)
                print(filename_plt + " created")
