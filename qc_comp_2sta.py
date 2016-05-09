#!/usr/bin/env python
#------------------------------------------------------------------------------
# Filename: qc_comp_2sta.py
#  Version: 1.0
#  Purpose: compare the quality between 2 stations/channels during a given period
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
from qc import *
from station_dictionnary import *
from instruments import *

# PARAMETERS
# ----------
# filename of the 2 pkl files
filename_pkl1 = "/Users/bonaime/Documents/info/test_site/PKL/G.OBP.00.BHN.pkl"
filename_pkl2 = "/Users/bonaime/Documents/info/test_site/PKL/G.OBP.10.BHN.pkl"

# Time period
period = (UTCDateTime(2015, 9, 10), UTCDateTime(2015, 12, 30))
# ----------

Hour = arange(0, 23, 1)
Hour = arange(0, 23)

try:
    fi1 = open(filename_pkl1, 'r')
    fi2 = open(filename_pkl2, 'r')
    S1 = cPickle.load(fi1)
    S2 = cPickle.load(fi2)
    if S1.sampling_rate != S2.sampling_rate:
        print "-----------------------------------------------------------------"
        print "!! CAREFUL !! sampling frequencies are different !!"
        print "!! interpolation of PSD frequency vectors not implemented yet  !!"
        print "!! Result show might be wrong !!"
        print "----------------------------------------------------------------"

except:
    print "No pickle file found for one of the 2 stations"
else:
    print "Use the pickle files : " + filename_pkl1 + " and " + filename_pkl2

    P = arange(period[0], period[1], 3600)
    TP = [UTCDateTime(t.year, t.month, t.day, t.hour) for t in P]
    T1H = [UTCDateTime(t.year, t.month, t.day, t.hour) for t in S1.times_used]
    T2H = [UTCDateTime(t.year, t.month, t.day, t.hour) for t in S2.times_used]

    TPs = intersect1d(TP, T1H)
    TPs = intersect1d(TPs, T2H)

    i1 = searchsorted(T1H, TPs)
    i2 = searchsorted(T2H, TPs)
    ip = searchsorted(TP, TPs)

    TP = TP[ip[0]:(ip[-1] + 1)]
    ip = ip - ip[0]
    TPnum = date2num(TP)

    coverage = empty((len(TP), 2))
    coverage.fill(nan)
    coverage[ip, :] = 1
    coverage = ma.array(coverage, mask=isnan(coverage))

    PSD1 = S1.psd[i1, :]
    PSD2 = S2.psd[i2, :]
    mpsd1 = median(PSD1, axis=0)
    mpsd2 = median(PSD2, axis=0)
    stdmpsd1 = std(PSD1, axis=0)
    stdmpsd2 = std(PSD2, axis=0)
    PSDdiff = empty((len(TP), len(S1.period_bins)))
    PSDdiff.fill(nan)
    for i, p in enumerate(S1.period_bins):
        PSDdiff[ip, i] = PSD2[:, i] - PSD1[:, i]
    PSDdiff = ma.array(PSDdiff, mask=isnan(PSDdiff))

    Spikes1 = S1.spikes[i1]
    Spikes2 = S2.spikes[i2]
    Spikesdiff = empty((len(TP)))
    Spikesdiff.fill(nan)
    Spikesdiff[ip] = Spikes2 - Spikes1
    Spikesdiff = ma.array(Spikesdiff, mask=isnan(Spikesdiff))

    HourUsed = array([t.hour for t in TPs])
    print HourUsed
    psdH1 = array([array(PSD1[HourUsed == h, :]).mean(axis=0) for h in Hour])
    spH1 = array([array(Spikes1[HourUsed == h, :]).mean(axis=0) for h in Hour])
    psdH2 = array([array(PSD2[HourUsed == h, :]).mean(axis=0) for h in Hour])
    spH2 = array([array(Spikes2[HourUsed == h, :]).mean(axis=0) for h in Hour])

    # Plot

    fig = figure(figsize=(8, 10), facecolor='w', edgecolor='k')
    ax_coverage = fig.add_axes([0.1, 0.08, 0.6, 0.03])
    ax_mpsd = fig.add_axes([0.15, 0.62, 0.7, 0.3])
    ax_diff_psd = fig.add_axes([0.1, 0.23, 0.6, 0.3])
    ax_diff_mpsd_hour = fig.add_axes([0.71, 0.23, 0.2, 0.3])
    ax_colorbar_diff_mpsd_hour = fig.add_axes([0.65, 0.54, 0.25, 0.015])
    ax_diff_spikes = fig.add_axes([0.1, 0.12, 0.6, 0.1])
    ax_diff_spikes_hour = fig.add_axes([0.71, 0.12, 0.2, 0.1])

    # Coverage
    cmap = cm.gray
    cmap.set_bad('w')
    tempo = array((0, 1))
    tx, y = meshgrid(TPnum, tempo)
    hcov = ax_coverage.pcolormesh(tx, y, coverage.T, cmap=cmap)
    hcov.set_clim((0, 2))
    ax_coverage.set_xlim((min(TPnum), max(TPnum)))
    setp(ax_coverage.get_yticklabels(), visible=False)
    ax_coverage.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax_coverage.xaxis.set_major_formatter(DateFormatter('%D'))
    for label in ax_coverage.get_xticklabels():
        label.set_ha("right")
        label.set_rotation(30)
    ax_coverage.xaxis.set_label_position('bottom')

    # PSD
    ax_mpsd.errorbar(S1.period_bins, mpsd1, yerr=stdmpsd1,
                     color='k', linewidth=2, ecolor='k', elinewidth=0.2)
    ax_mpsd.errorbar(S1.period_bins, mpsd2, yerr=stdmpsd2,
                     color='m', linewidth=2, ecolor='m', elinewidth=0.2)
    ax_mpsd.fill_between(S1.period_bins, mpsd1, y2=mpsd2,
                         where=(mpsd1 > mpsd2), color='b', alpha=0.5)
    ax_mpsd.fill_between(S1.period_bins, mpsd1, y2=mpsd2,
                         where=(mpsd1 < mpsd2), color='r', alpha=0.5)
    # Noise models
    model_periods, high_noise = get_NHNM()
    ax_mpsd.plot(model_periods, high_noise, '0.7', linewidth=2)
    model_periods, low_noise = get_NLNM()
    ax_mpsd.plot(model_periods, low_noise, '0.7', linewidth=2)
    ax_mpsd.semilogx()
    ax_mpsd.axis((1e-1, 3e2, -190, -90))
    ax_mpsd.grid('on')
    ax_mpsd.set_ylabel('median power (dB)')
    ax_mpsd.set_xlabel('period (s)')
    ax_mpsd.set_title(
        "station %s (magenta) - station %s (black)" % (S2.id, S1.id))

    # SPECTRO diff
    cmap = cm.jet
    cmap.set_bad('w')
    tx, F = meshgrid(TPnum, 1. / S1.per_octaves)
    hh = ax_diff_psd.pcolormesh(tx, F, PSDdiff.T, cmap=cmap)
    ax_diff_psd.semilogy()
    hh.set_clim((-10, 10))
    ax_diff_psd.set_ylim((3e-3, 10))
    ax_diff_psd.set_xlim((min(TPnum), max(TPnum)))
    ax_diff_psd.xaxis.set_ticks(ax_coverage.xaxis.get_ticklocs())
    ax_diff_psd.set_ylabel('frequency (Hz)')
    ax_diff_psd.xaxis.set_ticklabels("")
    ax_diff_psd.grid('on')

    # Spikes diff
    ax_diff_spikes.plot(TPnum[Spikesdiff > 0], Spikesdiff[
                        Spikesdiff > 0], color='r', marker='.', linestyle='None', alpha=0.5)
    ax_diff_spikes.plot(TPnum[Spikesdiff <= 0], Spikesdiff[
                        Spikesdiff <= 0], color='b', marker='.', linestyle='None', alpha=0.5)
    ax_diff_spikes.set_ylim((-10, 10))
    ax_diff_spikes.set_xlim((min(TPnum), max(TPnum)))
    ax_diff_spikes.xaxis.set_ticks(ax_coverage.xaxis.get_ticklocs())
    ax_diff_spikes.grid('on')
    ax_diff_spikes.set_ylabel('detection / hour')
    ax_diff_spikes.xaxis.set_ticklabels("")

    # Spectro diff / Hour
    H24, F = np.meshgrid(Hour, 1. / S1.per_octaves)
    h = ax_diff_mpsd_hour.pcolormesh(H24, F, (psdH2 - psdH1).T, cmap=cmap)
    h.set_clim(-10, 10)
    ax_diff_mpsd_hour.semilogy()
    ax_diff_mpsd_hour.axis((0, 23, 3e-3, 10))
    ax_diff_mpsd_hour.grid('on')
    ax_diff_mpsd_hour.set_xticks(arange(0, 23, 4))
    ax_diff_mpsd_hour.yaxis.set_ticks_position('right')
    ax_diff_mpsd_hour.yaxis.set_label_position('right')
    ax_diff_mpsd_hour.xaxis.set_ticklabels("")

    cb = colorbar(h, cax=ax_colorbar_diff_mpsd_hour,
                  orientation='horizontal', ticks=linspace(-8, 8, 5), format='%i')
    cb.set_label("dB")
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')
    cb.ax.xaxis.label.set_position((1, .2))
    cb.ax.xaxis.label.set_verticalalignment('top')

    # Spikes diff / Hour
    ax_diff_spikes_hour.plot(Hour, spH2 - spH1, color=(0.5, 0.5, 0.5))
    ax_diff_spikes_hour.fill_between(
        Hour, spH2 - spH1, y2=0, where=(spH2 < spH1), color='b', alpha=0.5)
    ax_diff_spikes_hour.fill_between(
        Hour, spH2 - spH1, y2=0, where=(spH2 > spH1), color='r', alpha=0.5)
    ax_diff_spikes_hour.axis((0, 23, -8, 8))
    ax_diff_spikes_hour.grid('on')
    ax_diff_spikes_hour.set_xlabel('Hour (UTC)')
    ax_diff_spikes_hour.set_xticks(arange(0, 23, 4))
    ax_diff_spikes_hour.set_yticks(arange(-6, 7, 3))
    ax_diff_spikes_hour.yaxis.set_ticks_position('right')
    ax_diff_spikes_hour.yaxis.set_label_position('right')

    show()


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
