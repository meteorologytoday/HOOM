import matplotlib as mplt
mplt.use('Agg')

import matplotlib.pyplot as plt
import Nio, sys, argparse
import numpy as np
from scipy import signal
from pprint import pprint

def SpectralVariance(y):
    c = np.fft.rfft(y, norm="ortho") # Ortho means normalized by sqrt N
    return abs(c)**2.0
 
parser = argparse.ArgumentParser()
parser.add_argument('--input-dir')
parser.add_argument('--output-dir')
parser.add_argument('--res')
parser.add_argument('--casenames')
parser.add_argument('--data-file')
parser.add_argument('--varname')
parser.add_argument('--normalize', default="yes")
parser.add_argument('--logy_max', type=float, default=3)
parser.add_argument('--logy_min', type=float, default=-3)
parser.add_argument('--ylabel', default="")

args = parser.parse_args()

pprint(args)

casenames = args.casenames.split(",")

print("Going to compare these models:")
pprint(casenames)

tss = []
sps = []

new_casenames = []

for i in range(len(casenames)):

    try:

        f = Nio.open_file("%s/%s_%s/%s" % (args.input_dir, args.res, casenames[i], args.data_file), "r")

    except Exception as e:
    
        print("Error happens when doing casename %s. Going to ignore this one..." % casenames[i])

        continue
    
    
    new_casenames.append(casenames[i])

    ts = f.variables[args.varname][:]

    if args.normalize == "yes":
        ts = signal.detrend(ts)
        ts /= np.std(ts)

    sp = SpectralVariance(ts)

    tss.append(ts)
    sps.append(sp)

    f.close()

casenames = new_casenames

N = len(tss[0])
freq = np.fft.rfftfreq(N, d=1.0/12.0)
time = np.arange(N) / 12
nyears = N / 12.0
period = 1.0 / freq

marked_periods   = np.array([0.5, 1, 2, 3, 4, 5, 6,7,8,9, 10, 15, 20, 30, 40])

fig, ax = plt.subplots(3, 1, figsize=(15, 10))

ax[0].set_title("%s (%d years)" % (args.varname, nyears, ))
ax[0].set_xlabel("Time [years]")
ax[0].set_ylabel(args.ylabel)
ax[0].grid()
for i in range(len(casenames)): 
    ax[0].plot(time, tss[i], linewidth=2, label=casenames[i])

ax[0].legend()

ax[1].set_title("Spectrum Analysis")
ax[1].set_xlabel("Period [years]")
ax[1].set_ylabel("Intensity  $| \\hat{c}(\\omega) |^2$")
ax[1].grid()

for i in range(len(casenames)): 
    ax[1].loglog(period[1:], sps[i][1:], linewidth=2, label=casenames[i])

ax[1].legend()

ax[1].set_xticks(marked_periods)
ax[1].set_xticklabels([("%.1f" if v<1 else "%d") % (v,) for v in marked_periods])

ax[1].set_ylim([10**args.logy_min, 10**args.logy_max])


ax[2].set_title("Spectrum Analysis")
ax[2].set_xlabel("Period [years]")
ax[2].set_ylabel("Intensity  $| \\hat{c}(\\omega) |^2$")
ax[2].grid()

for i in range(len(casenames)): 
    ax[2].plot(np.log(period[1:]), sps[i][1:], linewidth=2, label=casenames[i])

ax[2].legend()

ax[2].set_xticks(np.log(marked_periods))
ax[2].set_xticklabels([("%.1f" if v<1 else "%d") % (v,) for v in marked_periods])

#ax[2].set_ylim([10**args.logy_min, 10**args.logy_max])


fig.savefig("%s/%s_mc_climate_indices_%s.png" % (args.output_dir, args.res, args.varname), dpi=200)

#plt.show()
