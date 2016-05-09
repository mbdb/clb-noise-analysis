# Summary of Sensors and digitizers responses
# Approximate !!! To be used for quick deconvolution based on pz representation
# For sensors :
# Sensitivity in V/m/s (usualy call "gain" in RESP files)
# Poles and Zeros in rad/s for a velocity input in m/s
# Gain is k (A0) factor, sometimes manually determine as 1/max of the flat response
# For digitizers :
# lsb is least significant byte
# Note : Only the sensitivity is given. Do not depend on anti-alias filters or other manufacterer's dependant element.
# Information usually taken from the IRIS NRL (http://www.iris.edu/NRL), or from other sources (lithoscope, manufacturers webpages, ...).
# !!! NEED CONSTANT CHECK !!!
# -------------------------------------------------------------------
# SENSORS
# -------------------------------------------------------------------

# STREIKEISEN
# Based on a sensitivity of 1500 V/M/S
# sts2 : generic response (high pass filter). More complex representation
# depend on generation
sts2 = {'gain': 1.2982e8,
        'sensitivity': 1500,
        'poles': [(-3.70E-02 + 3.70E-02j), (-3.70E-02 - 3.70E-02j), (-1.0804E+02 + 4.198E+02j), (-1.0804E+02 - 4.198E+02j), (-6.91E+02 + 0j)],
        'zeros': [0j, 0j]}
# Specific for CEA stations (why more poles ???)
sts2_cea = {'gain': 1.2982e8,
            'sensitivity': 1500,
            'poles': [(-3.70E-02 + 3.70E-02j), (-3.70E-02 - 3.70E-02j), (-2.51E+02 + 0j), (-1.31E+02 + 4.67E+02j), (-1.31E+02 - 4.67E+02j)],
            'zeros': [0j, 0j]}
# sts1
sts1 = {'gain': 3.948580E+03,
        'sensitivity': 2400.,
        'poles': [(-1.234000e-02 + 1.234000e-02j), (-1.234000e-02 - 1.234000e-02j), (-3.918000e+01 + 4.912000e+01j), (-3.918000e+01 - 4.912000e+01j)],
        'zeros': [0j, 0j]}

# GURALP
# Based on the lowest sensitivity (800 V/M/S or 1500 V/M/S)
# cmg40t (30s) - from NRL
cmg40t30s = {'gain': 5.71508E+08,
             'sensitivity': 800,
             'poles': [(-1.486E-01 + 1.486E-01j), (-1.486E-01 - 1.486E-01j), (-5.0265E+02 + 0j), (-1.005E+03 + 0j), (-1.131E+03 + 0j)],
             'zeros': [0j, 0j]}
# cmg40t (30s) -- from guralp website
# cmg40t30s = {'gain':  2304000.,
#	'sensitivity' : 800,
#        'poles': [(-2.356E-02+2.356E-02j),(-2.356E-02-2.356E-02j),(-8.00E+01+0j),(-1.40E+02+0j),(-1.60E+02+0j)],
#        'zeros': [0j, 0j]}
# from NRL
cmg40t60s = {'gain': 5.71508E+08,
             'sensitivity': 800,
             'poles': [(-7.4016E-02 + 7.4016E-02j), (-7.4016E-02 - 7.4016E-02j), (-5.0265E+02 + 0j), (-1.005E+03 + 0j), (-1.131E+03 + 0j)],
             'zeros': [0j, 0j]}
# cmg3esp - PY46
cmg3py46 = {'gain': 2.30426E+06,
            'sensitivity': 2000,
            'poles': [(-7.61E-03 - 7.61E-03j), (-7.61E-03 + 7.61E-03j), (-8.0E+01 + 0j), (-1.6E+02 + 0j), (-1.8E+02 + 0j)],
            'zeros': [0j, 0j]}
# cmg40 default
cmg40t = cmg40 = cmg40t30s

# cmg3esp (30s) - similar to CMG40T except sensitivity -
cmg3esp30s = cmg40t30s
cmg3esp30s['sensitivity'] = 1500.
# cmg3esp (60s) - similar to CMG40T except sensitivity -
cmg3esp60s = cmg40t60s
cmg3esp60s['sensitivity'] = 2000.
# CMG3ESP (120s) - from NRL
cmg3esp120s = {'gain': 5.71508E+08,
               'sensitivity': 2000,
               'poles': [(-3.7008E-02 + 3.7008E-02j), (-3.7008E-02 - 3.7008E-02j), (-5.0265E+02 + 0j), (-1.005E+03 + 0j), (-1.131E+03 + 0j)],
               'zeros': [0j, 0j]}
# cmg3esp default
cmg3esp = cmg3esp120s

# cmg3esp - PY46
cmg3py46 = {'gain': 5.71508E+08,
            'sensitivity': 2000,
            'poles': [(-7.61E-03 + 7.61E-03j), (-7.61E-03 - 7.61E-03j), (-8.0E+01 + 0j), (-1.6E+02 + 0j), (-1.8E+02 + 0j)],
            'zeros': [0j, 0j]}

# cmg3t (30s) - similar to cmg3esp
cmg3t30s = cmg3esp30s
# cmg3t (60s) - similar to cmg3esp
cmg3t60s = cmg3esp60s
# cmg3t (120s) - similar to cmg3esp
cmg3t120s = cmg3esp120s
# cmg3t default
cmg3t = cmg3t120s

# AGECODAGIS
# noemax (20s) - From the Agecodagis website - values for the
# 034-0097-0.PZ - BUT DEPENDS ON the Sensor ID !!!!
noemax = {'gain': 1.,
          'sensitivity': 1500.,
          'poles': [(-0.2934739120327773 + 0j), (-0.2958623648278818 + 0j), (-3.9355748904513392 + 0j), (-216.96877577294066 + 0j)],
          'zeros': [0j, 0j, (-190.55241144076527 + 0j), (-3.6304229442730045 + 0j)]}

# KINEMETRICS
# ss-1
# From Helene
ss1 = {'gain': 1.,
       'sensitivity': 1.7066E+04,
       'poles': [(-2.22 + 2.22j), (-2.22 - 2.22j)],
       'zeros': [0j, 0j]}
# From NRL
# ss1 = {'gain': 1.,
#	'sensitivity' : 345.,
#        'poles': [(-4.44+4.44j), (-4.44-4.44j)],
#        'zeros': [0j, 0j]}

# LENNARTZ
# le3d1s
# From NRL
le3d1s = {'gain': 1.,
          'sensitivity': 400.,
          'poles': [(-4.4214 + 4.66j), (-4.4214 - 4.66j), (-2.105 + 0j)],
          'zeros': [0j, 0j, 0j]}
# Le3d5s
# From Lennartz technical documentation
le3d5s = {'gain': 1.,
          'sensitivity': 400,
          'poles': [(-0.888 + 0.888j), (-0.888 - 0.888j), (-0.220 + 0j)],
          'zeros': [0j, 0j, 0j]}
# From NRL
# le3d5s = {'gain': 1.,
#	'sensitivity' : 400,
#        'poles': [(-0.85+0.87j), (-0.85-0.87j), (-0.427+0j)],
#        'zeros': [0j, 0j, 0j]}
# Le3d20s
le3d20s = {'gain': 1.,
           'sensitivity': 1000.,
           'poles': [(-0.22 + 0.235j), (-0.22 - 0.235j), (-0.23 + 0j)],
           'zeros': [0j, 0j, 0j]}

# SERCEL
# From Helene - Do not correspond to NRL - CHECK
l4c3d = {'gain': 1.,
         'sensitivity': +1.15613E+04,
         'poles': [(-4.442 + 4.442j), (-4.442 - 4.442j)],
         'zeros': [0j, 0j]}

# NANOMETRICS
# Trillium 40
t40 = {'gain': 1.104000E+05,
       'sensitivity': 1553.,
       'poles': [(-2.41E+02 - 1.78E+02j), (-2.41E+02 + 1.78E+02j), (-5.35E+02 - 7.19E+02j), (-5.35E+02 + 7.19E+02j), (-8.63E+01 + 0j), (-1.103E-01 + 1.11E-01j), (-1.103E-01 - 1.11E-01j)],
       'zeros': [0j, 0j, -6.88E+01, -3.23E+02, -2.53E+03]}
# Trillium Compact
# From Iris/Passcal website; modified according to Maxime (15/5/11)
# tcompact = {'gain' : 7.4452e+11,
#	'sensitivity' : 749.1,
#	'poles' : [(-3.691E-02+3.712E-02j), (-3.691E-02-3.712E-02j), (-3.371E+02), (-3.739E+02+4.755E+02j), (-3.739E+02-4.755E+02j), (-5.884E+02+1.508E+03j), (-5.884E+02-1.508E+03j)],
#	'zeros' : [0j, 0j, -4.341E+02]}
# From NRL
tcompact = {'gain': 8.184000E+11,
            'sensitivity': 749.1,
            'poles': [(-3.691E-02 + 3.712E-02j), (-3.691E-02 - 3.712E-02j), (-3.371E+02), (-3.739E+02 + 4.755E+02j), (-3.739E+02 - 4.755E+02j), (-5.884E+02 + 1.508E+03j), (-5.884E+02 - 1.508E+03j)],
            'zeros': [0j, 0j, -4.341E+02]}
# Trillium 120
t120 = {'gain': 3.080E+05,
        'sensitivity': 1201.,
        'poles': [(-3.852E-02 + 3.658E-02j), (-3.852E-02 - 3.658E-02j), (-178. + 0j), (-135. + 160.j), (-135. - 160.j), (-671. + 1154.j), (-671. - 1154.j)],
        'zeros': [0j, 0j, -90., -160.7, -3108.]}
# Trillium 120PH (from NRL)
t120ph = {'gain': 8.318710E+17,
          'sensitivity': 1202.5,
          'poles': [(-3.6614e-02 + 3.7059e-02j), (-3.6614e-02 - 3.7059e-02j), (-3.255e+01), (-1.42e+02), (-3.64e+02 + 4.04e+02j), (-3.64e+02 - 4.04e+02j), (-1.26e+03), (-4.9e+03 + 5.204e+03j), (-4.9e+03 - 5.204e+03j), (-7.1e+03 + 1.7e+03j), (-7.1e+03 - 1.7e+03j)],
          'zeros': [0j, 0j, -31.63, -160., -350., -3177.]}
t120pa = t120
# Trillium 240 (!!! carefull second generation. Sensitivity is 1196.6 for
# first generation)
t240 = {'gain': 4.517000E+05,
        'sensitivity': 1.168200e+03,
        'poles': [(-1.770000e-02 + 1.760000e-02j), (-1.770000e-02 - 1.760000e-02j), (-1.267000e+02), (-1.920000e+02 + 2.591000e+02j), (-1.920000e+02 - 2.591000e+02j), (-5.577000e+02 + 1.143000e+03j), (-5.577000e+02 - 1.143000e+03j)],
        'zeros': [0j, 0j, -9.166000e+01, -1.601000e+02, -3.207000e+03]}

# -------------------------------------------------------------------
# DIGITIZERS
# -------------------------------------------------------------------

# Agecodagis !!!  TO CHECK !!! / not from NRL except k6
titan3xt = {'lsb': 745E-09}
titan3xtold = {'lsb': 298E-09}
titan3nt = {'lsb': 2.373E-06}
titan6nt = {'lsb': 2.373E-06}
titan_cs5321 = {'lsb': 839E-09}  # From RESP file for ANTF
titan_cs5323 = {'lsb': 596E-09}
titan_hi7190 = {'lsb': 729E-09}
k6bb = {'lsb': 9.31E-09}
k6 = {'lsb': 5.36000471680415e-07}
k3 = k6
kephren = k6
# Quanterra
q4120 = {'lsb': 2.3506944656659914e-06}
q330sr = {'lsb': 2.384188064754548e-06}
# !! Only for channels 1-3. channels 4-6 equivalent to q330sr
q330hr = {'lsb': 5.96047016188637e-07}
q330 = q330sr  # default for Q330
q330s = q330  # default for Q330S
# Nanometrics
# !!! 40V peak to peak correspond to a gain of 0.4
taurus40pp = {'lsb': 2.5000075000225e-06}
# !!! 16V peak to peak correspond to a gain of 1
taurus16pp = {'lsb': 1.000003000009e-06}
# !!! 8V peak to peak correspond to a gain of 2
taurus8pp = {'lsb': 5.000015000045e-07}
# !!! 4V peak to peak correspond to a gain of 4
taurus4pp = {'lsb': 2.5000075000225e-07}
taurus = taurus16pp  # default for taurus
trident40pp = {'lsb': 2.5000075000225e-06}  # Same as taurus
trident16pp = {'lsb': 1.000003000009e-06}
trident8pp = {'lsb': 5.000015000045e-07}
trident4pp = {'lsb': 2.5000075000225e-07}
trident = trident40pp  # default for trident
centaur = {'lsb': 2.5000075000225e-07}
# Reftek
reftek130a = {'lsb': 1.589499132928223e-06}
rt130a = reftek130a
# Others
tcheque = {'lsb': 1.2e-06}  # Digitizers used at Savonnieres from Plomerova
cea = {'lsb': 1.912e-06}
# geostar={'lsb' : 3.0518e-04/4.0816} # Ajout du pont diviseur de tension
geostar = {'lsb': 3.0518e-04}
acqstrf = {'lsb': 1.5259e-04}
