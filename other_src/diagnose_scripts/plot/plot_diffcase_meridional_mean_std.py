import matplotlib as mplt
mplt.use('Agg')

from matplotlib import rc

default_linewidth = 2.0;
default_ticksize = 10.0;

mplt.rcParams['lines.linewidth'] =   default_linewidth;
mplt.rcParams['axes.linewidth'] =    default_linewidth;
mplt.rcParams['xtick.major.size'] =  default_ticksize;
mplt.rcParams['xtick.major.width'] = default_linewidth;
mplt.rcParams['ytick.major.size'] =  default_ticksize;
mplt.rcParams['ytick.major.width'] = default_linewidth;

#rc('font', **{'family':'sans-serif', 'serif': 'Bitstream Vera Serif', 'sans-serif': 'MS Reference Sans Serif', 'size': 20.0});
rc('font', **{'size': 20.0});
rc('axes', **{'labelsize': 20.0});
rc('mathtext', **{'fontset':'stixsans'});
rc(('xtick.major','ytick.major'), pad=20)

#import matplotlib.font_manager as fm;
#print("%s: %d"%(fm.FontProperties().get_name(),fm.FontProperties().get_weight()));

import matplotlib.pyplot as plt

import sys, argparse
from netCDF4 import Dataset
import numpy as np
from pprint import pprint

parser = argparse.ArgumentParser()
parser.add_argument('--data-dir')
parser.add_argument('--ref-data-dir')
parser.add_argument('--data-file')

parser.add_argument('--casenames')
parser.add_argument('--ref-casenames')
parser.add_argument('--legends')

parser.add_argument('--domain-file')
parser.add_argument('--indexing', default=":")
parser.add_argument('--varname-mean')
parser.add_argument('--varname-var')
parser.add_argument('--display-varname')
parser.add_argument('--yscale', type=float, default=1.0)
parser.add_argument('--ylabel', default="")
parser.add_argument('--extra-title', default="")
parser.add_argument('--extra-filename', default="")
parser.add_argument('--colors')
parser.add_argument('--linestyles', type=str)
parser.add_argument('--y-offset', type=float, default=0.0)
parser.add_argument('--ymax-mean', type=float)
parser.add_argument('--ymax-std', type=float)
parser.add_argument('--tick-levs-mean', type=float)
parser.add_argument('--tick-levs-std', type=float)

args = parser.parse_args()

pprint(args)

casenames  = args.casenames.split(",")
legends    = args.legends.split(",")
colors     = args.colors.split(",")
linestyles = args.linestyles.split(",")
indexing   = args.indexing.split(",")
figsize   = tuple(map(float, args.figsize.split(',')))

indexes = []
print("Constructing indexing")
for i, content in enumerate(indexing):
    if content == ":":
        indexes.append(slice(None))
    else:
        indexes.append(int(content))

indexes = tuple(indexes)
print("Indices: ", indexes)    
    
print("Going to compare these models:")
pprint(casenames)

datas1 = []
datas2 = []

for i in range(N_cases):

    print("Loading case: %s" % casenames[i])    
    print("        case: %s" % ref_casenames[i])    
    f1 = Dataset("%s/%s/%s" % (args.data_dir, casenames[i], args.data_file), "r")
    f2 = Dataset("%s/%s/%s" % (args.ref_data_dir, ref_casenames[i], args.data_file), "r")
    
    datas1.append([f1.variables[args.varname_mean][idxes] / scale , f1.variables[args.varname_std][idxes] / scale])
    datas2.append([f2.variables[args.varname_mean][idxes] / scale, f2.variables[args.varname_std][idxes] / scale])

    f1.close()
    f2.close()


f = Dataset(args.domain_file, "r")
lat = f.variables["yc"][:, 1]
f.close()

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)

ax1.set_title("Mean of %s%s" % (args.display_varname, args.extra_title))
ax2.set_xlabel("Latitude [deg]")
ax2.set_xticks([-90, -60, -30, 0, 30, 60, 90])

ax1.set_ylabel(args.ylabel)
ax2.set_ylabel(args.ylabel)


for i in range(N_cases):
    _mean = datas1[i][0] - datas2[i][0]
    _std  = datas1[i][1] - datas2[i][1]
    ax1.plot(lat, _mean, linewidth=2, label=legends[i], color=colors[i], linestyle=linestyles[i])
    ax2.plot(lat, _std,  linewidth=2, label=legends[i], color=colors[i], linestyle=linestyles[i])

ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
fig.subplots_adjust(right=0.7, bottom=0.2)


ax1.grid(True)
ax2.grid(True)

ax1.set_ylim(np.array([-1.0, 1.0]) * args.ymax_mean_diff)
ax2.set_ylim(np.array([ 0.0, 1.0]) * args.ymax_std_diff)

fig.savefig("%s/diffcase_meridional_mean_std_%s_%s%s.png" % (args.output_dir, args.varname_mean, args.varname_var, args.extra_filename), dpi=200)

#plt.show()
