#!/usr/bin/env python
#------------------------------------------------------------------------------
# Filename: qc_comp_2periods.py
#  Version: 1.0
#  Purpose: compare the quality of a station/channel between 2 different periods of time
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
# Various paths
PATH_PKL = "/home/jerome/RESIF/CLB/Validation_Site/QC/PKL/XX/"
# PATH_PLT="/home/jerome/RESIF/CLB/Validation_Site/QC/Temp/"

# List of stations
STA = ['PY41']

CHAN = ['HHE']
# Time span
period1 = (UTCDateTime(2012, 7, 18), UTCDateTime(2012, 8, 1))
period2 = (UTCDateTime(2013, 1, 10), UTCDateTime(2013, 1, 23))

# ----------

Hour = arange(0, 23, 1)
# Loop over stations
for sta in STA:
    net = eval(sta)['network']
    locid = eval(sta)['locid']
    # Loop over channels
    for chan in CHAN:
        # Try to load a pickle file for this sta/net/chan
        filename_pkl = PATH_PKL + net + "." + sta + "." + locid + "." + chan + ".pkl"
        try:
            fi = open(filename_pkl, 'r')
            S = cPickle.load(fi)
        except:
            print "No pickle file found for " + net + "." + sta + "." + locid + "." + chan
        else:
            print "Use the pickle file : " + filename_pkl

            times_used = array(S.times_used)
            times_data = array(S.times_data)
            starttime1 = max(min(times_used), period1[0])
            endtime1 = min(max(times_used), period1[1])
            period1 = (starttime1, endtime1)
            starttime2 = max(min(times_used), period2[0])
            endtime2 = min(max(times_used), period2[1])
            period2 = (starttime2, endtime2)

            bool_times_select = (times_used > starttime1) & (
                times_used < endtime1)
            times_used1 = times_used[bool_times_select]
            times_data1 = times_data[
                (times_data[:, 0] > starttime1) & (times_data[:, 1] < endtime1)]
            psd1 = S.psd[bool_times_select, :]
            mpsd1 = median(psd1, axis=0)
            spikes1 = S.spikes[bool_times_select]

            bool_times_select = (times_used > starttime2) & (
                times_used < endtime2)
            times_used2 = times_used[bool_times_select]
            times_data2 = times_data[
                (times_data[:, 0] > starttime2) & (times_data[:, 1] < endtime2)]
            psd2 = S.psd[bool_times_select, :]
            mpsd2 = median(psd2, axis=0)
            spikes2 = S.spikes[bool_times_select]

            HourUsed1 = array([t.hour for t in times_used1])
            psdH1 = array([array(psd1[HourUsed1 == h, :]).mean(axis=0)
                           for h in Hour])
            spH1 = array([array(spikes1[HourUsed1 == h, :]).mean(axis=0)
                          for h in Hour])
            HourUsed2 = array([t.hour for t in times_used2])
            psdH2 = array([array(psd2[HourUsed2 == h, :]).mean(axis=0)
                           for h in Hour])
            spH2 = array([array(spikes2[HourUsed2 == h, :]).mean(axis=0)
                          for h in Hour])

            H24, P = np.meshgrid(Hour, S.per_octaves)

            # figure
            fig = figure(figsize=(8, 10), facecolor='w', edgecolor='k')
            ax_coverage_1 = fig.add_axes([0.13, 0.4, 0.02, 0.45])
            ax_coverage_2 = fig.add_axes([0.20, 0.4, 0.02, 0.45])
            ax_mpsd = fig.add_axes([0.3, 0.6, 0.6, 0.3])
            ax_diff_mpsd = fig.add_axes([0.3, 0.35, 0.6, 0.24])
            ax_diff_mpsd_hour = fig.add_axes([0.3, 0.1, 0.6, 0.24])
            ax_diff_spikes_hour = fig.add_axes([0.1, 0.1, 0.18, 0.24])
            ax_colorbar_diff_mpsd_hour = fig.add_axes([0.7, 0.05, 0.2, 0.02])

            # plots
            # coverage
            for start, end in times_data1:
                ax_coverage_1.axhspan(
                    start.datetime, end.datetime, color='k', lw=0)

            starts2 = date2num(times_used2)
            ends2 = date2num(times_used2 + S.len)
            for start, end in times_data2:
                ax_coverage_2.axhspan(
                    start.datetime, end.datetime, color='m', lw=0)

            ymin = min(min(times_data1[:, 0]), min(times_data2[:, 0]))
            ymax = max(max(times_data1[:, 1]), max(times_data2[:, 1]))
            ax_coverage_1.axis((0, 1, ymin.datetime, ymax.datetime))
            ax_coverage_2.axis((0, 1, ymin.datetime, ymax.datetime))
            setp(ax_coverage_1.get_xticklabels(), visible=False)
            setp(ax_coverage_2.get_xticklabels(), visible=False)
            setp(ax_coverage_2.get_yticklabels(), visible=False)
            ax_coverage_1.yaxis.set_major_locator(mdates.AutoDateLocator())
            ax_coverage_1.yaxis.set_major_formatter(DateFormatter('%D'))
            for label in ax_coverage_1.get_yticklabels():
                label.set_va("top")
                label.set_rotation(45)
            ax_coverage_1.set_ylabel("Time 1 (reference) = %5.1f days (black)" % (
                diff(times_data1).sum() / 86400))
            ax_coverage_2.set_ylabel("Time 2 = %5.1f days (magenta)" % (
                diff(times_data2).sum() / 86400))
            ax_coverage_1.yaxis.set_label_position('right')
            ax_coverage_2.yaxis.set_label_position('right')

            # PSD
            ax_mpsd.fill_between(S.period_bins, mpsd1, y2=mpsd2, where=(
                mpsd1 > mpsd2), color='b', alpha=0.5)
            ax_mpsd.fill_between(S.period_bins, mpsd1, y2=mpsd2, where=(
                mpsd1 < mpsd2), color='r', alpha=0.5)
            ax_mpsd.plot(S.period_bins, mpsd1, 'k', linewidth=2)
            ax_mpsd.plot(S.period_bins, mpsd2, 'm', linewidth=2)

            # Noise models
            model_periods, high_noise = get_NHNM()
            ax_mpsd.plot(model_periods, high_noise, '0.7', linewidth=2)
            model_periods, low_noise = get_NLNM()
            ax_mpsd.plot(model_periods, low_noise, '0.7', linewidth=2)

            ax_mpsd.semilogx()
            ax_mpsd.axis((1e-1, 5e2, -190, -90))
            ax_mpsd.grid('on')
            ax_mpsd.set_ylabel('median power (dB)')
            ax_mpsd.xaxis.set_ticks_position('top')
            ax_mpsd.xaxis.set_label_position('top')
            ax_mpsd.yaxis.set_ticks_position('right')
            ax_mpsd.yaxis.set_label_position('right')

            ax_mpsd.set_title(S.id)
            ax_mpsd.title.set_y(1.1)

            ax_diff_mpsd.fill_between(
                S.period_bins, mpsd2 - mpsd1, y2=0, where=(mpsd2 > mpsd1), color='r', alpha=0.5)
            ax_diff_mpsd.fill_between(
                S.period_bins, mpsd2 - mpsd1, y2=0, where=(mpsd2 < mpsd1), color='b', alpha=0.5)
            ax_diff_mpsd.semilogx()
            ax_diff_mpsd.grid('on')
            ax_diff_mpsd.axis((1e-1, 5e2, -12, 12))
            setp(ax_diff_mpsd.get_xticklabels(), visible=False)
            ax_diff_mpsd.set_xlabel('Periods (s)')
            ax_diff_mpsd.set_ylabel('T2 - T1 (dB)')
            ax_diff_mpsd.yaxis.set_ticks_position('right')
            ax_diff_mpsd.yaxis.set_label_position('right')

            h = ax_diff_mpsd_hour.pcolor(P, H24, (psdH2 - psdH1).T)
            h.set_clim(-10, 10)
            ax_diff_mpsd_hour.semilogx()
            ax_diff_mpsd_hour.axis((1e-1, 5e2, 0, 23))
            ax_diff_mpsd_hour.grid('on')
            ax_diff_mpsd_hour.set_ylabel('Hour (UTC)')
            ax_diff_mpsd_hour.set_xlabel('Periods (s)')
            cb = colorbar(h, cax=ax_colorbar_diff_mpsd_hour,
                          orientation='horizontal', ticks=linspace(-8, 8, 5), format='%i')
            cb.set_label("T2 -T1 (dB)")
            ax_diff_mpsd_hour.set_yticks(arange(0, 23, 4))
            ax_diff_mpsd_hour.yaxis.set_ticks_position('right')
            ax_diff_mpsd_hour.yaxis.set_label_position('right')

            ax_diff_spikes_hour.plot(spH2 - spH1, Hour, color=(0.5, 0.5, 0.5))
            ax_diff_spikes_hour.fill_betweenx(
                Hour, spH2 - spH1, x2=0, where=(spH2 < spH1), color='b', alpha=0.5)
            ax_diff_spikes_hour.fill_betweenx(
                Hour, spH2 - spH1, x2=0, where=(spH2 > spH1), color='r', alpha=0.5)
            ax_diff_spikes_hour.axis((-8, 8, 0, 23))
            ax_diff_spikes_hour.grid('on')
            ax_diff_spikes_hour.set_ylabel('Hour (UTC)')
            ax_diff_spikes_hour.set_xlabel('delta_detections')
            ax_diff_spikes_hour.set_yticks(arange(0, 23, 4))
            ax_diff_spikes_hour.set_xticks(arange(-6, 7, 3))

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
