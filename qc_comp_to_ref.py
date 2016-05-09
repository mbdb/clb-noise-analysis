#!/usr/bin/env python
#------------------------------------------------------------------------------
# Filename: qc_comp_to_ref.py
#  Version: 1.0
#  Purpose: compare the median PSD at a bunch of stations with NLNM/NHNM and A/B sites reference curves
#     Note: QC pkl files are needed
#   Author: Jerome Vergne
#    Email: jerome.vergne@unsitra.fr
#
#------------------------------------------------------------------------------
import os
import glob
import cPickle
from obspy.core import *
from obspy.xseed import Parser
from scipy import interpolate
from qc import *
from station_dictionnary import *

# PARAMETERS
# ----------
# Various paths
PATH_PKL = "/home/jerome/RESIF/CLB/Validation_Site/QC/PKL/XX/"
PATH_PLT_OUT = "/home/jerome/RESIF/CLB/Validation_Site/QC/PLT/XX/"
NOISE_MODEL_FILE = "/home/jerome/RESIF/CLB/Validation_Site/QC/data/noise_models.npz"
CLASS_MODEL_FILE = "/home/jerome/RESIF/CLB/Validation_Site/QC/data/class_models.npz"

# List of stations
# stations permanentes
#STA=['TRBF'  , 'ANTF' , 'ARBF' , 'ATE' , 'CALF' ,'CFF' , 'CHIF' ,'CHMF' ,'DOU' ,'FNEB' ,'ISO' , 'LRVF' , 'MLS' , 'MONQ' , 'OG02' ,  'OG35' ,'OGAG' ,'OGDI' ,'OGSM' , 'PYLO' , 'RENF' , 'RSL' ,  'RUSF' , 'SAOF' ,'SJAF' ,'SURF' , 'ECH' ,'SSB' , 'CLF' , 'LOT2' ,'MTT2' ,'ORT2' ,'PGT2' ,'ROT2']
# sites de test
#STA=['NAN1' , 'AGO' , 'AJAC' , 'OG03' ,  'OGMO' , 'OSSF' , 'PLDF' , 'VIEF' , 'PY48' , 'PY48_3' , 'PY48_4' , 'ARZ' , 'AUDR2' , 'AUDR' , 'BARI2' , 'BARI3' , 'BARI' , 'CAVD' , 'CLER2' , 'CLER' , 'CLEV' , 'CRMS' , 'CRUG' , 'GPD' , 'GRAU3' , 'GRAUV' , 'GRN1' , 'GRN5' , 'MYY' , 'N25E' , 'N25I' , 'NEUF' , 'POLI' , 'POUR' , 'PY41' , 'ROST' , 'ROVI' , 'SACE' , 'SAND' , 'SAVO' , 'STMU'  , 'URBES' , 'VAUD' , 'VERTU' , 'VILS' , 'VINC' , 'WING']
# quelques sites fermes
STA = ['MYY', 'POLI', 'ROVI', 'CAVD', 'URBES', 'BARI2', 'BARI3',  'VAUD',
       'CRUG', 'VERTU', 'GPD', 'SAVO', 'SACE', 'GRN5', 'OGMO', 'ROST', 'NEUF']
# Quelques sites ouverts
#STA=['VILS' , 'AGO' , 'SAND' , 'PLDF' , 'CRMS' , 'N25I' , 'GRN1' , 'ARZ' , 'OG03' , 'CLEV' , 'N25E' , 'OSSF' , 'VIEF' , 'FILF' , 'CLER' , 'CLER2' , 'AUDR' , 'AUDR2' , 'GRAUV' , 'NAN1' , 'WING']

COMP = ['Z', 'N', 'E']

# Time span for the estimation of the median PSD
period = (UTCDateTime(2012, 1, 1), UTCDateTime(2012, 12, 31))

# Period ranges used to sort stations
# One png file per period range (and per component) will be generated
Per_sort = ([50., 200.], [.2, .5], [5., 20.])

# Plotting arguments
show_class = True  # To compare to A/B sites theoretical curves
show_txt_class = True  # indicate if a station is A/B on the right of the plot
# dB above class limit to accept a site as A or B (useful if PSd median is
# just about the class limit in a restricted period range
allow_above_class = 3.
Clim = ((0, 50))  # colorscale limits compared to NLNM

basename_output = "sort_CLOS"  # basename for output png graphs

# ---- END PARAMETERS ----
# ------------------------

# should be modified acording to Per_sort and freq sampling
Per = logspace(log10(0.1), log10(200), 100)

y = arange(0, len(STA) + 1)

per_ref, psd_ref = get_NLNM()

aa = per_ref.argsort()
per_ref = per_ref[aa]
psd_ref = psd_ref[aa]
finterp = interpolate.interp1d(per_ref, psd_ref, bounds_error=False)
psd_ref = finterp(Per)

cmap = cm.jet
cmap.set_bad('w')

if show_class:
    PSDclass = empty([3, len(Per)])
    per_refA, psd_refA, per_refB, psd_refB = get_class()
    finterp = interpolate.interp1d(per_refA, psd_refA, bounds_error=False)
    PSDclass[0, :] = finterp(Per) - psd_ref
    finterp = interpolate.interp1d(per_refB, psd_refB, bounds_error=False)
    PSDclass[1, :] = finterp(Per) - psd_ref
    PSDclass[2, :] = finterp(Per) - psd_ref
    PSDclass = ma.array(PSDclass, mask=isnan(PSDclass))
    yclass = arange(0, 3)


for comp in COMP:
    PSD = empty([len(STA), len(Per)])
    PSD.fill(nan)

    Times_used = [UTCDateTime(t.year, t.month, t.day, t.hour)
                  for t in arange(period[0], period[1], 3600)]
    for ista, sta in enumerate(STA):
        net = eval(sta)['network']
        locid = eval(sta)['locid']
        # Try to load a pickle file for this sta/net/chan
        filename_pkl = PATH_PKL + net + "." + sta + "." + locid + ".??" + comp + ".pkl"
        # filename_pkl=PATH_PKL+"*."+sta+".*."+comp+".pkl"
        try:
            fi = open(glob.glob(filename_pkl)[0], 'r')
            S = cPickle.load(fi)
        except:
            print "No pickle file found for " + filename_pkl
            sys.stdout.flush()
        else:
            print "Read the pickle file : " + glob.glob(filename_pkl)[0]
            sys.stdout.flush()

            try:
                periods, psd = S.get_percentile(
                    percentile=50, hist_cum=None, time_lim=period)
            except:
                print S.id + " have no data in the time range"
            else:
                finterp = interpolate.interp1d(
                    periods, psd, bounds_error=False)
                PSD[ista, :] = finterp(Per) - psd_ref
                #SPIKES[ista]=mean([S.spikes[i] for i,t in enumerate(S.times_used) if(t>period[0] and t<period[1])])

    PSD = ma.array(PSD, mask=isnan(PSD))

    for per_sort in Per_sort:
        print per_sort
        iper_sort = searchsorted(Per, per_sort)
        isort = (PSD[:, iper_sort[0]:iper_sort[-1]].mean(axis=1)).argsort()
        PSDs = PSD[isort, :]
        PSDs = append(PSDs, [PSDs[-1, :]], axis=0)
        PSDs = ma.array(PSDs, mask=isnan(PSDs))
        # isort=isort[-1::-1]
        ysta = array(STA)
        ysta = ysta[isort]
        ysta = append(ysta, ysta[-1])
        # Figure
        fig = figure(figsize=(8, 10), facecolor='w', edgecolor='k')
        if show_class:
            ax_psd = fig.add_axes([0.1, 0.2, 0.8, 0.75])
            ax_colorbar = fig.add_axes([0.65, 0.14, 0.25, 0.015])
            ax_class = fig.add_axes([0.1, 0.04, 0.8, 0.05])
        else:
            ax_psd = fig.add_axes([0.1, 0.15, 0.8, 0.8])
            ax_colorbar = fig.add_axes([0.65, 0.08, 0.25, 0.015])

        h = ax_psd.pcolormesh(Per, y, PSDs, cmap=cmap)
        ax_psd.semilogx()
        h.set_clim(Clim)
        ax_psd.set_yticks(y + 0.5)
        ax_psd.set_ylim((0, y[-1] + 1))
        ax_psd.yaxis.set_ticklabels(ysta)
        hl = ax_psd.plot([per_sort[0], per_sort[0]], [y[0], y[-1]], 'm-')
        hl = ax_psd.plot([per_sort[-1], per_sort[-1]], [y[0], y[-1]], 'm-')
        ax_psd.set_xlim((Per[0], Per[-1]))
        ax_psd.grid('on')
        ax_psd.set_xlabel('period [s]')
        ax_psd.set_title(
            "PSD - NLNM / Comp %s / sorted by T=[%4.1f - %4.1f] s" % (comp, per_sort[0], per_sort[-1]))
        cb = colorbar(h, cax=ax_colorbar, orientation='horizontal',
                      ticks=arange(Clim[0], Clim[1] + 1, 10), format='%i')
        cb.set_label("dB above NLNM")

        if show_class:
            hclass = ax_class.pcolormesh(Per, yclass, PSDclass, cmap=cmap)
            ax_class.semilogx()
            hclass.set_clim(Clim)
            ax_class.set_yticks(yclass + 0.5)
            ax_class.set_ylim((0, yclass[-1]))
            ax_class.yaxis.set_ticklabels(['class A', 'class B'])
            ax_class.set_xlim((Per[0], Per[-1]))
            ax_class.grid('on')
            if show_txt_class:
                for i, aa in enumerate(PSDs[:, iper_sort[0]:iper_sort[-1]]):
                    isA = np.all(
                        aa < (PSDclass[0, iper_sort[0]:iper_sort[-1]] + allow_above_class))
                    if isA:
                        txtclass = " A / "
                    else:
                        txtclass = " - / "
                    isB = np.all(aa < PSDclass[1, iper_sort[0]:iper_sort[-1]])
                    if isB:
                        txtclass += "B"
                    else:
                        txtclass += "-"
                    if i == (shape(PSDs)[0] - 1):
                        txtclass = ""

                    htxtclass = ax_psd.text(Per[-1], i + 0.5, txtclass)
                    htxtclass.set_verticalalignment('center')

        savefig(PATH_PLT_OUT + basename_output + "_%s_%4.1f_%4.1f.png" %
                (comp, per_sort[0], per_sort[-1]))


def get_NLNM():
    """
    Returns periods and psd values for the New Low Noise Model.
    """
    data = np.load(NOISE_MODEL_FILE)
    periods = data['model_periods']
    nlnm = data['low_noise']
    return (periods, nlnm)


def get_NHNM():
    """
    Returns periods and psd values for the New High Noise Model.
    """
    data = np.load(NOISE_MODEL_FILE)
    periods = data['model_periods']
    nhnm = data['high_noise']
    return (periods, nhnm)


def get_class():
    """
    Returns periods and psd values for the New High Noise Model.
    """
    data = np.load(CLASS_MODEL_FILE)
    periodsA = data['perA']
    nmA = data['dbA']
    periodsB = data['perB']
    nmB = data['dbB']
    return (periodsA, nmA, periodsB, nmB)
