#!/usr/bin/env python
#------------------------------------------------------------------------------
# Filename: qc.py
#  Version: 1.0
#  Purpose: 	Clone of psd.py from obspy to estimate the quality of a site
#		by computing the spectrogram, PPSD and number of detections per hour
#   Author: Jerome Vergne, Bonaime Sebastien
#    Email: jerome.vergne@unsitra.fr, bonaime@ipgp.fr
#
#------------------------------------------------------------------------------
# -*- coding: utf-8 -*-

"""
Various Routines Related to Spectral Estimation
"""
from __future__ import with_statement
import os
import warnings
import cPickle
import math
import bisect
import argparse
import glob
from sys import stdout, exit
import numpy as np
from pylab import arange, array, DateFormatter, transpose, colorbar, linspace,\
setp, zeros, size, ma, cm, NaN, vstack, hstack
import matplotlib.dates as mdates
#from obspy.core import *
from obspy.core import Trace, Stream, UTCDateTime
from obspy.clients.filesystem.sds import Client
#from obspy.core.util import MATPLOTLIB_VERSION
from obspy.signal.filter import bandpass
from obspy.signal.trigger import recursive_sta_lta, trigger_onset
from obspy.signal.invsim import cosine_taper
from obspy.signal.util import prev_pow_2
from obspy.io.xseed import Parser
from obspy.signal.spectral_estimation import get_nhnm, get_nlnm
from obspy.imaging.cm import pqlx, get_cmap

from instruments import *
from station_dictionnary import *
from IPython import embed

MATPLOTLIB_VERSION = "Exist"

if MATPLOTLIB_VERSION is None:
    # if matplotlib is not present be silent about it and only raise the
    # ImportError if matplotlib actually is used (currently in psd() and
    # PPSD())
    msg_matplotlib_ImportError = "Failed to import matplotlib. While this " \
        "is no dependency of obspy.signal it is however necessary for a " \
        "few routines. Please install matplotlib in order to be able " \
        "to use e.g. psd() or PPSD()."
    # set up two dummy functions. this makes it possible to make the docstring
    # of psd() look like it should with two functions as default values for
    # kwargs although matplotlib might not be present and the routines
    # therefore not usable

    def detrend_none():
        pass

    def window_hanning():
        pass
else:
    # Import matplotlib routines. These are no official dependency of
    # obspy.signal so an import error should really only be raised if any
    # routine is used which relies on matplotlib (at the moment: psd, PPSD).
    from matplotlib import mlab
    import matplotlib.pyplot as plt
    from matplotlib.dates import date2num
    from matplotlib.ticker import FormatStrFormatter
    from matplotlib.mlab import detrend_none, window_hanning

# PSD PARAMTERS
# -------------

CLASS_MODEL_FILE = os.path.join(os.path.dirname(__file__),
                                "data", "class_models.npz")
# do not change these variables, otherwise results may differ from PQLX!
PPSD_LENGTH = 3600  # psds are calculated on 1h long segments
# PPSD_STRIDE = 1800 # psds are calculated overlapping, moving 0.5h ahead
PPSD_STRIDE = 3600  # psds are calculated without overlapping
T1 = UTCDateTime(1970, 1, 1)
T2 = UTCDateTime(2030, 1, 1)

# DETECTION PARAMTERS (DP)
# ------------------------
# This has to be confirmed based on usual values used in observatories
DP_f_min = 2.
DP_f_max = 8.
DP_f_corners = 4
DP_sta = 1.
DP_lta = 5.
DP_max_len = 3.
DP_thres_high = 2.5
DP_thres_low = 2.0

# CMAPS
spectro_cmap = 'jet'


def psd(x, NFFT=256, Fs=2, detrend=detrend_none, window=window_hanning, noverlap=0):
    """
    Wrapper for `matplotlib.mlab.psd`_.
    """
    # check if matplotlib is available, no official dependency for obspy.signal
    if MATPLOTLIB_VERSION is None:
        raise ImportError(msg_matplotlib_ImportError)

    # check matplotlib version
    elif MATPLOTLIB_VERSION >= [0, 98, 4]:
        new_matplotlib = True
    else:
        new_matplotlib = False
    # build up kwargs that do not change with version 0.98.4
    kwargs = {}
    kwargs['NFFT'] = NFFT
    kwargs['Fs'] = Fs
    kwargs['detrend'] = detrend
    kwargs['window'] = window
    kwargs['noverlap'] = noverlap
    # add additional kwargs to control behavior for matplotlib versions higher
    # than 0.98.4. These settings make sure that the scaling is already done
    # during the following psd call for newer matplotlib versions.
    if new_matplotlib:
        kwargs['pad_to'] = None
        kwargs['sides'] = 'onesided'
        kwargs['scale_by_freq'] = True
    # do the actual call to mlab.psd
    Pxx, freqs = mlab.psd(x, **kwargs)
    # do scaling manually for old matplotlib versions
    if not new_matplotlib:
        Pxx = Pxx / Fs
        Pxx[1:-1] = Pxx[1:-1] * 2.0
    return Pxx, freqs


def fft_taper(data):
    """
    Cosine taper, 10 percent at each end.
    """
    data *= cosine_taper(len(data), 0.2)
    return data


def welch_taper(data):
    """
    Applies a welch window to data
    """
    data *= welch_window(len(data))
    return data


def welch_window(N):
    """
    Return a welch window for data of length N
    """
    n = math.ceil(N / 2.0)
    taper_left = np.arange(n, dtype=np.float64)
    taper_left = 1 - np.power(taper_left / n, 2)
    # first/last sample is zero by definition
    if N % 2 == 0:
        # even number of samples: two ones in the middle, perfectly symmetric
        taper_right = taper_left
    else:
        # odd number of samples: still two ones in the middle, however, not
        # perfectly symmetric anymore. right side is shorter by one sample
        nn = n - 1
        taper_right = np.arange(nn, dtype=np.float64)
        taper_right = 1 - np.power(taper_right / nn, 2)
    taper_left = taper_left[::-1]
    # first/last sample is zero by definition
    taper_left[0] = 0.0
    taper_right[-1] = 0.0
    taper = np.concatenate((taper_left, taper_right))
    return taper

# ------------------------------------------------------------
# ------------------------------------------------------------


class QC(object):
    """
    Class to compile probabilistic power spectral densities+STA/LTA detections
    for one combination of network/station/location/channel/sampling_rate.
    """

    def __init__(self, stats, paz=None, dataless=None, skip_on_gaps=False):
        """
        Initialize the PPSD object setting all fixed information on the station
        that should not change afterwards to guarantee consistent spectral
        estimates.
        The necessary instrument response information can be provided in two
        ways:

        * Providing a dataless file. This is the safer way but it might a bit
          slower because for every processed time segment the response
          information is extracted from the parser.
        * Providing a dictionary containing poles and zeros information. Be
          aware that this leads to wrong results if the instrument's response
          is changing with data added to the PPSD. Use with caution!

        :type stats: :class:`~obspy.core.trace.Stats`
        :param stats: Stats of the station/instrument to process
        :type paz: dict (optional)
        :param paz: Response information of instrument. If not specified the
                information is supposed to be present as stats.paz.
        :type dataless: String (optional)
        :param dataless: Dataless file  with response information
        :type skip_on_gaps: Boolean (optional)
        :param skip_on_gaps: Determines whether time segments with gaps should
                be skipped entirely. McNamara & Buland merge gappy
                traces by filling with zeros. This results in a clearly
                identifiable outlier psd line in the PPSD visualization. Select
                `skip_on_gaps=True` for not filling gaps with zeros which might
                result in some data segments shorter than 1 hour not used in
                the PPSD.
        """
        # check if matplotlib is available, no official dependency for
        # obspy.signal
        if MATPLOTLIB_VERSION is None:
            raise ImportError(msg_matplotlib_ImportError)

        if paz is not None and dataless is not None:
            msg = "Both paz and parser specified. Using parser object for " \
                  "metadata."
            warnings.warn(msg)

        self.id = "%(network)s.%(station)s.%(location)s.%(channel)s" % stats
        self.network = stats.network
        self.station = stats.station
        self.location = stats.location
        self.channel = stats.channel
        self.sampling_rate = stats.sampling_rate
        self.delta = 1.0 / self.sampling_rate
        # trace length for one hour piece
        self.len = int(self.sampling_rate * PPSD_LENGTH)
        # set paz either from kwarg or try to get it from stats
        self.paz = paz
        self.dataless = dataless
        self.parser = Parser(dataless)
        if skip_on_gaps:
            self.merge_method = -1
        else:
            self.merge_method = 0
        # nfft is determined mimicing the fft setup in McNamara&Buland paper:
        # (they take 13 segments overlapping 75% and truncate to next lower
        #  power of 2)
        #  - take number of points of whole ppsd segment (currently 1 hour)
        self.nfft = PPSD_LENGTH * self.sampling_rate
        #  - make 13 single segments overlapping by 75%
        #    (1 full segment length + 25% * 12 full segment lengths)
        self.nfft = self.nfft / 4.0
        #  - go to next smaller power of 2 for nfft
        self.nfft = prev_pow_2(self.nfft)

        #  - use 75% overlap (we end up with a little more than 13 segments..)
        self.nlap = int(0.75 * self.nfft)
        self.times_used = []
        self.times = self.times_used
        self.times_data = []
        self.times_gaps = []
        self.hist_stack = None
        self.psd = []
        self.spikes = []
        self.__setup_bins()

    def __setup_bins(self):
        """
        Makes an initial dummy psd and thus sets up the bins and all the rest.
        Should be able to do it without a dummy psd..
        """
        dummy = np.ones(self.len)
        spec, freq = mlab.psd(
            dummy, self.nfft, self.sampling_rate, noverlap=self.nlap)

        # leave out first entry (offset)
        freq = freq[1:]

        per = 1.0 / freq[::-1]
        self.freq = freq
        self.per = per
        # calculate left/rigth edge of first period bin, width of bin is one
        # octave
        per_left = per[0] / 2
        per_right = 2 * per_left
        # calculate center period of first period bin
        per_center = math.sqrt(per_left * per_right)
        # calculate mean of all spectral values in the first bin
        per_octaves_left = [per_left]
        per_octaves_right = [per_right]
        per_octaves = [per_center]
        # we move through the period range at 1/8 octave steps
        factor_eighth_octave = 2**(1. / 8)
        # do this for the whole period range and append the values to our lists
        while per_right < per[-1]:
            per_left *= factor_eighth_octave
            per_right = 2 * per_left
            per_center = math.sqrt(per_left * per_right)
            per_octaves_left.append(per_left)
            per_octaves_right.append(per_right)
            per_octaves.append(per_center)
        self.per_octaves_left = np.array(per_octaves_left)
        self.per_octaves_right = np.array(per_octaves_right)
        self.per_octaves = np.array(per_octaves)

        self.period_bins = per_octaves
        # mid-points of all the period bins
        self.period_bin_centers = np.mean((self.period_bins[:-1],
                                           self.period_bins[1:]), axis=0)
        # set up the binning for the db scale
        self.spec_bins = np.linspace(-200, -50, 301, endpoint=True)

    def __sanity_check(self, trace):
        """
        Checks if trace is compatible for use in the current PPSD instance.
        Returns True if trace can be used or False if not.
        """
        if trace.id != self.id:
            return False
        if trace.stats.sampling_rate != self.sampling_rate:
            return False
        return True

    def __insert_used_time(self, utcdatetime):
        """
        Inserts the given UTCDateTime at the right position in the list keeping
        the order intact.
        Add sorting of self.psd and self.spikes arrays
        """
        if len(self.times_used) > 1:
            i = hstack((array(self.times_used), utcdatetime)).argsort()
            self.psd = self.psd[i]
            self.spikes = self.spikes[i]
        bisect.insort(self.times_used, utcdatetime)

    def __insert_gap_times(self, stream):
        """
        Gets gap information of stream and adds the encountered gaps to the gap
        list of the PPSD instance.
        """
        self.times_gaps += [[gap[4], gap[5]] for gap in stream.get_gaps()]

    def __insert_data_times(self, stream):
        """
        Gets gap information of stream and adds the encountered gaps to the gap
        list of the PPSD instance.
        """
        self.times_data += [[tr.stats.starttime, tr.stats.endtime]
                            for tr in stream]

    def __check_time_present(self, utcdatetime):
        """
        Checks if the given UTCDateTime is already part of the current PPSD
        instance. That is, checks if from utcdatetime to utcdatetime plus 1
        hour there is already data in the PPSD.
        Returns True if adding an one hour piece starting at the given time
        would result in an overlap of the ppsd data base, False if it is OK to
        insert this piece of data.
        """
        index1 = bisect.bisect_left(self.times_used, utcdatetime)
        index2 = bisect.bisect_right(self.times_used, utcdatetime + PPSD_LENGTH)
        if index1 != index2:
            return True
        return False

    def add(self, stream, verbose=True):
        """
        Process all traces with compatible information and add their spectral
        estimates to the histogram containg the probabilistic psd.
        Also ensures that no piece of data is inserted twice.
        """
        # return later if any changes were applied to the ppsd statistics
        changed = False
        # prepare the list of traces to go through
        if isinstance(stream, Trace):
            stream = Stream([stream])
        # select appropriate traces
        stream = stream.select(id=self.id,
                               sampling_rate=self.sampling_rate)
        # save information on available data and gaps
        self.__insert_data_times(stream)
        self.__insert_gap_times(stream)
        # merge depending on skip_on_gaps set during __init__
        stream.merge(self.merge_method, fill_value=0)

        for tr in stream:
            # the following check should not be necessary due to the select()..
            if not self.__sanity_check(tr):
                msg = "Skipping incompatible trace."
                warnings.warn(msg)
                continue
            t1 = tr.stats.starttime
            t2 = tr.stats.endtime
            while t1 + PPSD_LENGTH <= t2:
                if self.__check_time_present(t1):
                    msg = "Already computed time spans detected (e.g. %s), " + \
                          "skipping these slices."
                    msg = msg % t1
                    print msg
                else:
                    # throw warnings if trace length is different than one
                    # hour..!?!
                    slice = tr.slice(t1, t1 + PPSD_LENGTH)
                    success = self.__process(slice)
                    if success:
                        self.__insert_used_time(t1)
                        if verbose:
                            stdout.write("\r adding %s" % t1)
                            stdout.flush()
                        changed = True
                t1 += PPSD_STRIDE  # advance half an hour
        if verbose:
            stdout.write("\r")
            stdout.flush()
        return changed

    def __process(self, tr):
        """
        Processes a one-hour segment of data and adds the information to the
        PPSD histogram. If Trace is compatible (station, channel, ...) has to
        checked beforehand.
        """
        # XXX DIRTY HACK!!
        if len(tr) == self.len + 1:
            tr.data = tr.data[:-1]
        # one last check..
        if len(tr) != self.len:
            msg = "Got an non-one-hour piece of data to process. Skipping"
            warnings.warn(msg)
            print len(tr), self.len
            return False
        # being paranoid, only necessary if in-place operations would follow
        tr.data = tr.data.astype("float64")
        # if trace has a masked array we fill in zeros
        try:
            tr.data[tr.data.mask] = 0.0
        # if its no masked array, we get an AttributeError and have nothing to
        # do
        except AttributeError:
            pass

        # get instrument response preferably from parser object
        try:
            paz = self.parser.get_paz(self.id, datetime=tr.stats.starttime)
        except Exception, e:
            if self.dataless is not None:
                msg = "Error getting response from dataless:\n%s: %s\n" \
                      "Skipping time segment(s)."
                msg = msg % (e.__class__.__name__, e.message)
                warnings.warn(msg)
# -------------- CHANGED (comment next line)-----------
#                return False
            paz = self.paz
        if paz is None:
            msg = "Missing poles and zeros information for response " \
                  "removal. Skipping time segment(s)."
            warnings.warn(msg)
            return False

        tr.simulate(paz_remove=paz, remove_sensitivity=True,
                    paz_simulate=None, simulate_sensitivity=False)

        # Spike detector
        # band-pass filtering
        tr_filt = bandpass(tr.data, DP_f_min, DP_f_max,
                           self.sampling_rate, corners=DP_f_corners)
        # STA/LTA
        #cft = recStaltaPy(tr_filt,DP_sta*self.sampling_rate,DP_lta*self.sampling_rate)
        cft = recursive_sta_lta(tr_filt, int(
            DP_sta * self.sampling_rate), int(DP_lta * self.sampling_rate))

        on_of = trigger_onset(cft, DP_thres_high, DP_thres_low,
                              max_len=DP_max_len * self.sampling_rate)

        # go to acceleration
        tr.data = np.gradient(tr.data, self.delta)

        # use our own wrapper for mlab.psd to have consistent results on all
        # matplotlib versions
        spec, freq = psd(tr.data, self.nfft, self.sampling_rate,
                         detrend=mlab.detrend_linear, window=fft_taper,
                         noverlap=self.nlap)

        # leave out first entry (offset)
        spec = spec[1:]

        # working with the periods not frequencies later so reverse spectrum
        spec = spec[::-1]

        # go to dB
        spec = np.log10(spec)
        spec *= 10

        spec_octaves = []
        # do this for the whole period range and append the values to our lists
        for per_left, per_right in zip(self.per_octaves_left, self.per_octaves_right):
            spec_center = spec[(per_left <= self.per) &
                               (self.per <= per_right)].mean()
            spec_octaves.append(spec_center)
        spec_octaves = np.array(spec_octaves)

        hist, self.xedges, self.yedges = np.histogram2d(
            self.per_octaves, spec_octaves, bins=(self.period_bins, self.spec_bins))

        try:
            self.hist_stack += hist
            self.psd = vstack((self.psd, spec_octaves))
            self.spikes = hstack((self.spikes, array(len(on_of))))
        except TypeError:
            # initialize (only during first run)
            self.hist_stack = hist
            self.psd = spec_octaves
            self.spikes = array(len(on_of))
        return True

    def get_percentile(self, percentile=50, hist_cum=None, time_lim=(T1, T2)):
        """
        Returns periods and approximate psd values for given percentile value.

        :type percentile: int
        :param percentile: percentile for which to return approximate psd
                value. (e.g. a value of 50 is equal to the median.)
        :type hist_cum: `numpy.ndarray` (optional)
        :param hist_cum: if it was already computed beforehand, the normalized
                cumulative histogram can be provided here (to avoid computing
                it again), otherwise it is computed from the currently stored
                histogram.
        :returns: (periods, percentile_values)
        """
        if hist_cum is None:
            hist_cum = self.__get_normalized_cumulative_histogram(
                time_lim=time_lim)
        # go to percent
        percentile = percentile / 100.0
        if percentile == 0:
            # only for this special case we have to search from the other side
            # (otherwise we always get index 0 in .searchsorted())
            side = "right"
        else:
            side = "left"
        percentile_values = [col.searchsorted(percentile, side=side)
                             for col in hist_cum]
        # map to power db values
        percentile_values = self.spec_bins[percentile_values]
        return (self.period_bin_centers, percentile_values)

    def __get_normalized_cumulative_histogram(self, time_lim=(T1, T2)):
        """
        Returns the current histogram in a cumulative version normalized per
        period column, i.e. going from 0 to 1 from low to high psd values for
        every period column.
        """
        # sum up the columns to cumulative entries
        hist_cum = (self.__get_ppsd(time_lim=(T1, T2))).cumsum(axis=1)
        # normalize every column with its overall number of entries
        # (can vary from the number of self.times because of values outside
        #  the histogram db ranges)
        norm = hist_cum[:, -1].copy()
        # avoid zero division
        norm[norm == 0] = 1
        hist_cum = (hist_cum.T / norm).T
        return hist_cum

    def __get_psd(self, time_lim=(T1, T2), per_lim=(.1, 1)):
        """
        Returns the mean psd in the per_lim and time_lim range
        """
        bool_select_time = np.all([array(self.times_used) > time_lim[0],
                                   array(self.times_used) < time_lim[1]], axis=0)
        bool_select_per = np.all(
            [self.per_octaves > per_lim[0], self.per_octaves < per_lim[1]], axis=0)
        mpsd = self.psd[bool_select_time, :]
        mpsd = mpsd[:, bool_select_per]
        return mpsd.mean(1)

    def __get_ppsd(self, time_lim=(T1, T2)):
        """
        Returns the normalised ppsd in the time_lim range
        """
        bool_select_time = np.all([array(self.times_used) > time_lim[0],
                                   array(self.times_used) < time_lim[1]], axis=0)

        for i, l_psd in enumerate(self.psd[bool_select_time, :]):
            hist, xedges, yedges = np.histogram2d(
                self.per_octaves, l_psd, bins=(self.period_bins, self.spec_bins))
            try:
                hist_stack += hist
            except:
                hist_stack = hist
        hist_stack = hist_stack * 100.0 / i
        return hist_stack

    def save(self, filename):
        """
        Saves PPSD instance as a pickled file that can be loaded again using
        pickle.load(filename).
        """
        # with open(filename, "w") as file:
        #    cPickle.dump(self, file)

        cPickle.dump(self, open(filename, "wb"))


# -------------------- P L O T -----------------------------------------
    def plot(self, cmap, filename=None,
             starttime=T1, endtime=T2,
             show_percentiles=False, percentiles=[10, 50, 90],
             show_class_models=True, grid=True, title_comment=False):
        """
        Plot the QC resume figure
        If a filename is specified the plot is saved to this file, otherwise
        a plot window is shown.

        :type filename: str (optional)
        :param filename: Name of output file
        :type show_percentiles: bool (optional)
        :param show_percentiles: Enable/disable plotting of approximated
                percentiles. These are calculated from the binned histogram and
                are not the exact percentiles.
        :type percentiles: list of ints
        :param percentiles: percentiles to show if plotting of percentiles is
                selected.
        :type show_class_models: bool (optional)
        :param show_class_models: Enable/disable plotting of class models.
        :type grid: bool (optional)
        :param grid: Enable/disable grid in histogram plot.
        :type cmap: cmap
        :param cmap: Colormap for PPSD.
        """

        # COMMON PARAMETERS
        psd_db_limits = (-180, -110)
        psdh_db_limits = (-200, -90)
        f_limits = (5e-3, 20)
        per_left = (10, 1, .1)
        per_right = (100, 10, 1)
        # -----------------

        # Select Time window
        # -----------
        times_used = array(self.times_used)
        starttime = max(min(times_used), starttime)
        endtime = min(max(times_used), endtime)
        bool_times_select = (times_used > starttime) & (times_used < endtime)
        times_used = times_used[bool_times_select]
        psd = self.psd[bool_times_select, :]
        spikes = self.spikes[bool_times_select]

        hist_stack = self._QC__get_ppsd(time_lim=(starttime, endtime))

        Hour = arange(0, 23, 1)
        HourUsed = array([t.hour for t in times_used])
        Day_span = (endtime - starttime) / 86400.
        # -----------

        # FIGURE and AXES
        fig = plt.figure(figsize=(9.62, 13.60), facecolor='w', edgecolor='k')
        ax_ppsd = fig.add_axes([0.1, 0.68, 0.9, 0.28])
        ax_coverage = fig.add_axes([0.1, 0.56, 0.64, 0.04])
        ax_spectrogram = fig.add_axes([0.1, 0.31, 0.64, 0.24])
        ax_spectrogramhour = fig.add_axes([0.76, 0.31, 0.20, 0.24])
        ax_freqpsd = fig.add_axes([0.1, 0.18, 0.64, 0.12])
        ax_freqpsdhour = fig.add_axes([0.76, 0.18, 0.20, 0.12])
        ax_spikes = fig.add_axes([0.1, 0.05, 0.64, 0.12])
        ax_spikeshour = fig.add_axes([0.76, 0.05, 0.20, 0.12])

        ax_col_spectrogram = fig.add_axes([0.76, 0.588, 0.20, 0.014])
        ax_col_spectrogramhour = fig.add_axes([0.76, 0.57, 0.20, 0.014])

        ########################### COVERAGE
        ax_coverage.xaxis_date()
        ax_coverage.set_yticks([])
        # plot data coverage
        starts = date2num([a.datetime for a in times_used])
        ends = date2num([a.datetime for a in times_used + PPSD_LENGTH])
        for start, end in zip(starts, ends):
            ax_coverage.axvspan(start, end, 0, 0.7, alpha=0.5, lw=0)
        # plot data really available
        aa = [(start, end) for start, end in self.times_data if (
            (end - start) > PPSD_LENGTH)]  # avoid very small gaps	otherwise very long to plot
        for start, end in aa:
            start = date2num(start.datetime)
            end = date2num(end.datetime)
            ax_coverage.axvspan(start, end, 0.7, 1, facecolor="g", lw=0)
        # plot gaps
        aa = [(start, end) for start, end in self.times_gaps if (
            (end - start) > PPSD_LENGTH)]  # avoid very small gaps	otherwise very long to plot
        for start, end in aa:
            start = date2num(start.datetime)
            end = date2num(end.datetime)
            ax_coverage.axvspan(start, end, 0.7, 1, facecolor="r", lw=0)
        # Compute uncovered periods
        starts_uncov = ends[:-1]
        ends_uncov = starts[1:]
        # Keep only major uncovered periods
        ga = (ends_uncov - starts_uncov) > (PPSD_LENGTH) / 86400
        starts_uncov = starts_uncov[ga]
        ends_uncov = ends_uncov[ga]

        ax_coverage.set_xlim(starttime.datetime, endtime.datetime)

        # labels
        ax_coverage.xaxis.set_ticks_position('top')
        ax_coverage.tick_params(direction='out')

        ax_coverage.xaxis.set_major_locator(mdates.AutoDateLocator())
        if Day_span > 5:
            ax_coverage.xaxis.set_major_formatter(DateFormatter('%D'))
        else:
            ax_coverage.xaxis.set_major_formatter(DateFormatter('%D-%Hh'))
            for label in ax_coverage.get_xticklabels():
                label.set_fontsize(10)


        for label in ax_coverage.get_xticklabels():
            label.set_ha("right")
            label.set_rotation(-25)

        ########################### SPECTROGRAM
        ax_spectrogram.xaxis_date()
        t = date2num([a.datetime for a in times_used])
        f = 1. / self.per_octaves
        T, F = np.meshgrid(t, f)
        spectro = ax_spectrogram.pcolormesh(
            T, F, transpose(psd), cmap=spectro_cmap)
        spectro.set_clim(*psd_db_limits)

        spectrogram_colorbar = colorbar(spectro, cax=ax_col_spectrogram,
                                        orientation='horizontal',
                                        ticks=linspace(psd_db_limits[0],
                                                       psd_db_limits[1], 5),
                                        format='%i')
        spectrogram_colorbar.set_label("dB")
        spectrogram_colorbar.set_clim(*psd_db_limits)
        spectrogram_colorbar.ax.xaxis.set_ticks_position('top')
        spectrogram_colorbar.ax.xaxis.label.set_position((1.1, .2))
        spectrogram_colorbar.ax.yaxis.label.set_horizontalalignment('left')
        spectrogram_colorbar.ax.yaxis.label.set_verticalalignment('bottom')

        ax_spectrogram.grid(which="major")
        ax_spectrogram.semilogy()
        ax_spectrogram.set_ylim(f_limits)
        ax_spectrogram.set_xlim(starttime.datetime, endtime.datetime)
        ax_spectrogram.set_xticks(ax_coverage.get_xticks())
        setp(ax_spectrogram.get_xticklabels(), visible=False)
        ax_spectrogram.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        ax_spectrogram.set_ylabel('Frequency [Hz]')
        ax_spectrogram.yaxis.set_label_coords(-0.08, 0.5)

        ########################### SPECTROGRAM PER HOUR
        #psdH=array([array(psd[HourUsed==h,:]).mean(axis=0) for h in Hour])
        psdH = zeros((size(Hour), size(self.per_octaves)))
        for i, h in enumerate(Hour):
            a = array(psd[HourUsed == h, :])
            A = ma.masked_array(
                a, mask=~((a > psdh_db_limits[0]) & (a < psdh_db_limits[1])))
            psdH[i, :] = ma.getdata(A.mean(axis=0))
        psdH = array([psdH[:, i] - psdH[:, i].mean()
                      for i in arange(0, psdH.shape[1])])
        H24, F = np.meshgrid(Hour, f)
        spectroh = ax_spectrogramhour.pcolormesh(H24, F, psdH, cmap=cm.RdBu_r)
        spectroh.set_clim(-8, 8)

        spectrogram_per_hour_colorbar = colorbar(spectroh,
                                                 cax=ax_col_spectrogramhour,
                                                 orientation='horizontal',
                                                 ticks=linspace(-8, 8, 5),
                                                 format='%i')
        spectrogram_per_hour_colorbar.set_clim(-8, 8)

        ax_spectrogramhour.semilogy()
        ax_spectrogramhour.set_xlim((0, 23))
        ax_spectrogramhour.set_ylim(f_limits)
        ax_spectrogramhour.set_xticks(arange(0, 23, 4))
        ax_spectrogramhour.set_xticklabels(arange(0, 23, 4), visible=False)
        ax_spectrogramhour.yaxis.set_ticks_position('right')
        ax_spectrogramhour.yaxis.set_label_position('right')
        ax_spectrogramhour.yaxis.grid(True)
        ax_spectrogramhour.xaxis.grid(False)

        ########################### PSD BY PERIOD RANGE
        t = date2num([a.datetime for a in times_used])
        ax_freqpsd.xaxis_date()
        for pp in zip(per_left, per_right):
            mpsd = self._QC__get_psd(time_lim=(starttime, endtime), per_lim=pp)
            mpsdH = zeros(size(Hour)) + NaN
            for i, h in enumerate(Hour):
                a = array(mpsd[HourUsed == h])
                A = ma.masked_array(
                    a, mask=~((a > psdh_db_limits[0]) & (a < psdh_db_limits[1])))
                mpsdH[i] = ma.getdata(A.mean())
            ax_freqpsd.plot(t, mpsd)
            ax_freqpsdhour.plot(Hour, mpsdH - mpsdH.mean())
        ax_freqpsd.set_ylim(psd_db_limits)
        ax_freqpsd.set_xlim(starttime.datetime, endtime.datetime)
        ax_freqpsd.set_xticks(ax_coverage.get_xticks())
        setp(ax_freqpsd.get_xticklabels(), visible=False)
        ax_freqpsd.set_ylabel('Amplitude [dB]')
        ax_freqpsd.yaxis.set_label_coords(-0.08, 0.5)
        ax_freqpsd.yaxis.grid(False)
        ax_freqpsd.xaxis.grid(True)

        ########################### PSD BY PERIOD RANGE PER HOUR
        ax_freqpsdhour.set_xlim((0, 23))
        ax_freqpsdhour.set_ylim((-8, 8))
        ax_freqpsdhour.set_yticks(arange(-6, 7, 2))
        ax_freqpsdhour.set_xticks(arange(0, 23, 4))
        ax_freqpsdhour.set_xticklabels(arange(0, 23, 4), visible=False)
        ax_freqpsdhour.yaxis.set_ticks_position('right')
        ax_freqpsdhour.yaxis.set_label_position('right')

        ########################### SPIKES
        ax_spikes.xaxis_date()
        ax_spikes.bar(t, spikes, width=1. / 24)
        ax_spikes.set_ylim((0, 50))
        ax_spikes.set_xlim(starttime.datetime, endtime.datetime)
        ax_spikes.set_yticks(arange(10, 45, 10))
        ax_spikes.set_xticks(ax_coverage.get_xticks())
        #setp(ax_spikes.get_xticklabels(), visible=False)
        ax_spikes.set_ylabel("Detections [#/hour]")
        ax_spikes.yaxis.set_label_coords(-0.08, 0.5)
        ax_spikes.yaxis.grid(False)
        ax_spikes.xaxis.grid(True)

        # labels
        ax_spikes.xaxis.set_ticks_position('bottom')
        ax_spikes.tick_params(direction='out')
        ax_spikes.xaxis.set_major_locator(mdates.AutoDateLocator())
        if Day_span > 5:
            ax_spikes.xaxis.set_major_formatter(DateFormatter('%D'))
        else:
            ax_spikes.xaxis.set_major_formatter(DateFormatter('%D-%Hh'))
            for label in ax_spikes.get_xticklabels():
                label.set_fontsize(10)

        for label in ax_spikes.get_xticklabels():
            label.set_ha("right")
            label.set_rotation(25)

        ########################### SPIKES PER HOUR
        mspikesH = array([array(spikes[[HourUsed == h]]).mean() for h in Hour])
        ax_spikeshour.bar(Hour, mspikesH - mspikesH.mean(), width=1.)
        ax_spikeshour.set_xlim((0, 23))
        ax_spikeshour.set_ylim((-8, 8))
        ax_spikeshour.set_xticks(arange(0, 23, 4))
        ax_spikeshour.set_yticks(arange(-6, 7, 2))
        ax_spikeshour.set_ylabel("Daily variation")
        ax_spikeshour.set_xlabel("Hour [UTC]")
        ax_spikeshour.yaxis.set_ticks_position('right')
        ax_spikeshour.yaxis.set_label_position('right')
        ax_spikeshour.yaxis.set_label_coords(1.3, 1)

        ########################### plot gaps
        for start, end in zip(starts_uncov, ends_uncov):
            ax_spectrogram.axvspan(
                start, end, 0, 1, facecolor="w", lw=0, zorder=100)
            ax_freqpsd.axvspan(
                start, end, 0, 1, facecolor="w", lw=0, zorder=100)
            ax_spikes.axvspan(start, end, 0, 1,
                              facecolor="w", lw=0, zorder=100)

        # LEGEND
        leg = [str(xx) + '-' + str(yy) + ' s' for xx,
               yy in zip(per_left, per_right)]
        hleg = ax_freqpsd.legend(
            leg, loc=3, bbox_to_anchor=(-0.015, 0.75), ncol=size(leg))
        for txt in hleg.get_texts():
            txt.set_fontsize(8)

        # PPSD
        X, Y = np.meshgrid(self.xedges, self.yedges)
        ppsd = ax_ppsd.pcolormesh(X, Y, hist_stack.T, cmap=cmap)
        ppsd_colorbar = plt.colorbar(ppsd, ax=ax_ppsd)
        ppsd_colorbar.set_label("PPSD [%]")
        color_limits = (0, 30)
        ppsd.set_clim(*color_limits)
        ppsd_colorbar.set_clim(*color_limits)
        ax_ppsd.grid(b=grid, which="major")

        if show_percentiles:
            hist_cum = self.__get_normalized_cumulative_histogram(
                time_lim=(starttime, endtime))
            # for every period look up the approximate place of the percentiles
            for percentile in percentiles:
                periods, percentile_values = self.get_percentile(
                    percentile=percentile, hist_cum=hist_cum, time_lim=(starttime, endtime))
                ax_ppsd.plot(periods, percentile_values, color="black")

        # Noise models
        model_periods, high_noise = get_nhnm()
        ax_ppsd.plot(model_periods, high_noise, '0.4', linewidth=2)
        model_periods, low_noise = get_nlnm()
        ax_ppsd.plot(model_periods, low_noise, '0.4', linewidth=2)
        if show_class_models:
            classA_periods, classA_noise, classB_periods, classB_noise = get_class()
            ax_ppsd.plot(classA_periods, classA_noise, 'r--', linewidth=3)
            ax_ppsd.plot(classB_periods, classB_noise, 'g--', linewidth=3)

        ax_ppsd.semilogx()
        ax_ppsd.set_xlim(1. / f_limits[1], 1. / f_limits[0])
        ax_ppsd.set_ylim((-200, -80))
        ax_ppsd.set_xlabel('Period [s]')
        ax_ppsd.get_xaxis().set_label_coords(0.5, -0.05)
        ax_ppsd.set_ylabel('Amplitude [dB]')
        ax_ppsd.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))

        # TITLE
        title = "%s   %s -- %s  (%i segments)"
        title = title % (self.id, starttime.date, endtime.date,
                         len(times_used))
        if title_comment:
            fig.text(0.82, 0.978, title_comment, bbox=dict(
                facecolor='red', alpha=0.5), fontsize=15)

        ax_ppsd.set_title(title)
        # a=str(UTCDateTime().format_iris_web_service())
        plt.draw()

        if filename is not None:
            plt.savefig(filename)
            plt.close()
        else:
            plt.show()


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


def main():

    # Arguments
    argu_parser = argparse.ArgumentParser(
        description="Run the main code to compute quality sheets for a set of STATIONS")
    argu_parser.add_argument("-s", "--stations", nargs='*', required=True,
                             help='list of stations name Separate with spaces that must be in the station_dictionnary file ex: MEUD00 OBP10')
    argu_parser.add_argument("-b", "--starttime", default=UTCDateTime(2000, 1, 1), type=UTCDateTime,
                             help="Start time for processing. Various format accepted. Example : 2012,2,1 / 2012-02-01 / 2012,032 / 2012032 / etc ... See UTCDateTime for a complete list. Default is 2010-1-1")
    argu_parser.add_argument("-e", "--endtime", default=UTCDateTime(2055, 9, 16), type=UTCDateTime,
                             help="End time for processing. Various format accepted. Example : 2012,2,1 / 2012-02-01 / 2012,032 / 2012032 / etc ... See UTCDateTime for a complete list. Default is 2015-1-1")
    argu_parser.add_argument("-c", "--channels", nargs='+',
                             help="Process only CHANNELS. Do not use this option if you want to process all available channels indicated in station_dictionnary. Separate CHANNELS with spaces. No wildcard. Default is all channels")
    argu_parser.add_argument("-pkl", "--path_pkl", default='./PKL',
                             help="output directory for pkl files. Default is ./PKL  ")
    argu_parser.add_argument("-plt", "--path_plt", default='./PLT',
                             help="output directory for plt files. Default is ./PLT ")
    argu_parser.add_argument("-force_paz", default=False, action='store_true',
                             help="Use this option if you want don't want to use the dataless file specified in station_dictionnary. Only PAZ response computed from sensor and digitizer will be used. Use for debug only. Dataless recommended")
    argu_parser.add_argument("--color_map", default='pqlx',
                             help="Color map for PPSD. Default is pqlx")

    args = argu_parser.parse_args()

    # List of stations
    STA = args.stations
    chan_proc = args.channels

    # Time span
    start = args.starttime
    stop = args.endtime

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

    # Output Paths
    PATH_PKL = os.path.abspath(args.path_pkl)
    PATH_PLT = os.path.abspath(args.path_plt)

    # Create PKL and PLT directory if they don't exist

    if not os.path.exists(PATH_PKL):
        os.makedirs(PATH_PKL)

    if not os.path.exists(PATH_PLT):
        os.makedirs(PATH_PLT)

    # ----------

    # Check if all stations are in station_dictionnary

    for dict_station_name in STA:
        try:
            dict_station_name
        except:
            print dict_station_name + " is not in station_dictionnary.py"
            exit()

        try:
            eval(dict_station_name)['network']
        except:
            print dict_station_name + " does not have a network in station_dictionnary.py"
            exit()

        try:
            eval(dict_station_name)['station']
        except:
            print dict_station_name + " does not have a station in station_dictionnary.py"
            exit()
        try:
            eval(dict_station_name)['locid']
        except:
            print dict_station_name + " does not have a locid in station_dictionnary.py"
            exit()

        try:
            eval(dict_station_name)['dataless_file']
            # if there is a dataless_file, check if the file exists
            if not os.path.isfile(eval(dict_station_name)['dataless_file']):
                print dict_station_name + " does not have a valid dataless_file in station_dictionnary.py :"
                print eval(dict_station_name)['dataless_file'] + ' does not exist'
                exit()
        except:
            pass

        try:
            eval(dict_station_name)['path_data']
        except:
            print dict_station_name + " does not have a path_data in station_dictionnary.py"
            exit()

        if not os.path.exists(eval(dict_station_name)['path_data']):
            print dict_station_name + " does not have a valid path_data in station_dictionnary.py :"
            print eval(dict_station_name)['path_data'] + ' does not exist'
            exit()

    # Loop over stations
    for dict_station_name in STA:
        net = eval(dict_station_name)['network']
        sta = eval(dict_station_name)['station']
        locid = eval(dict_station_name)['locid']
        sds_path = eval(dict_station_name)['path_data']

        try:
            title_comment = eval(dict_station_name)['title_comment']
        except:
            title_comment = False
            pass

        # use only channels in station_dictionnay
        if isinstance(chan_proc, list):
            chan_proc = np.intersect1d(
                chan_proc, eval(dict_station_name)['channels'])
        else:
            chan_proc = eval(dict_station_name)['channels']

        # Look for a DATALESS
        dataless = None
        if not args.force_paz:
            try:
                dataless = glob.glob(eval(dict_station_name)['dataless_file'])[0]
            except:
                print "No valid dataless found for " + net + "." + sta + "." + locid
            else:
                print "Using dataless file : " + dataless
        paz = None
        if dataless is None:
            # Look for a PAZ
            try:
                sismo = eval((eval(dict_station_name)['sensor'].lower()))
                acq = eval((eval(dict_station_name)['digitizer'].lower()))

            except:
                print "No PAZ found for " + net + "." + sta + "." + locid
            else:
                paz = {'gain': sismo['gain'],
                       'poles': sismo['poles'],
                       'zeros': sismo['zeros'],
                       'sensitivity': sismo['sensitivity'] / acq['lsb']}
                print "PAZ from instruments.py: "
                print paz

        # exit if no dataless nor paz
        if dataless is None and paz is None:
            print "you must provide a dataless file or a sensor and a digitizer from instruments.py"
            exit()

        # exit if dataless AND paz are provided
        if dataless is not None and paz is not None:
            print "you must provide a dataless file or a sensor and a digitizer from instruments.py but not both !"
            exit()

        print "SDS archive is " + str(sds_path)

        # Loop over channels
        for chan in chan_proc:

            # Try to load a pickle file for this sta/net/chan
            filename_pkl = PATH_PKL + '/' + net + "." + \
                sta + "." + locid + "." + chan + ".pkl"
            is_pickle = False
            try:
                pkl_file = open(filename_pkl, 'r')
                quality_check = cPickle.load(pkl_file)
            except:
                print "No pickle file found for " + net + "." + sta + "." + locid + "." + chan + " (looked for " + filename_pkl + ")"
            else:
                print "Use the pickle file : " + filename_pkl
                is_pickle = True

            sds_client = Client(sds_path)

            print "Reading %s.%s.%s.%s in SDS archive from %s to %s" % (net, sta, locid, chan, start, stop)
            all_streams = sds_client.get_waveforms(
                net, sta, locid, chan, start, stop)

            if len(all_streams.traces) == 0:
                print "No data found for %s.%s.%s.%s in SDS archive from %s to %s" % (net, sta, locid, chan, start, stop)
                exit()

            print "Processing data"
            # Initiate the QC
            if is_pickle is False:
                quality_check = QC(all_streams[0].stats, dataless=dataless,
                                   paz=paz, skip_on_gaps=True)
                is_pickle = True

            # Remove the minutes before the next hour (useful when
            # computing statistics per hour)
            # mst = min start time
            mst = min([temp.stats.starttime for temp in all_streams])
            mst = UTCDateTime(mst.year, mst.month,
                              mst.day, mst.hour, 0, 0) + 3600.
            all_streams.trim(starttime=mst, nearest_sample=False)

            # Trim to stop when asked
            all_streams.trim(endtime=stop)

            # Remove the minutes after the last hour (useful when
            # computing statistics per hour)
            # met = max end time
            met = max([temp.stats.endtime for temp in all_streams])
            met = UTCDateTime(met.year, met.month, met.day, met.hour, 0, 0)
            all_streams.trim(endtime=met, nearest_sample=False)

            # Add all streams to QC
            quality_check.add(all_streams)

            # save
            # This can be 2 levels below to save a little bit of time
            # (save and loading)
            if is_pickle:
                print filename_pkl
                quality_check.save(filename_pkl)
                print filename_pkl + " updated"
            else:
                print "!!!!! Nothing saved/created for  " + net + "." + sta + "." + locid
            # Plot
            if is_pickle and len(quality_check.times_used) > 0:
                filename_plt = PATH_PLT + '/' + net + "." + \
                    sta + "." + locid + "." + chan + ".png"
                quality_check.plot(cmap=cmap, filename=filename_plt,
                                   show_percentiles=True, starttime=start,
                                   endtime=stop, title_comment=title_comment)
                print filename_plt + " updated"


if __name__ == '__main__':
    main()
