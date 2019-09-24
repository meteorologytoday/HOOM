import matplotlib as mplt
mplt.use('Agg')

import matplotlib.pyplot as plt
import sys, argparse
from netCDF4 import Dataset
import numpy as np
from pprint import pprint

parser = argparse.ArgumentParser()
parser.add_argument('--input-dir')
parser.add_argument('--output-dir')
parser.add_argument('--casenames')
parser.add_argument('--legends')
parser.add_argument('--data-file')
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
parser.add_argument('--y-range-mean', type=str, default="")
parser.add_argument('--y-range-std', type=str, default="")

args = parser.parse_args()

pprint(args)

casenames  = args.casenames.split(",")
legends    = args.legends.split(",")
colors     = args.colors.split(",")
linestyles = args.linestyles.split(",")
indexing   = args.indexing.split(",")

indices = []
print("Constructing indexing")
for i, content in enumerate(indexing):
    if content == ":":
        indices.append(slice(None))
    else:
        indices.append(int(content))

indices = tuple(indices)
print("Indices: ", indices)    
    
if args.y_range_mean != "":
    args.y_range_mean = [float(n) for n in args.y_range_mean.split(",")]

if args.y_range_std != "":
    args.y_range_std = [float(n) for n in args.y_range_std.split(",")]



print("Going to compare these models:")
pprint(casenames)

new_casenames = []
mean_datas = []
std_datas = []
for i in range(len(casenames)):

    try:

        f = Dataset("%s/%s/%s" % (args.input_dir, casenames[i], args.data_file), "r")

    except Exception as e:
   
        print(str(e)) 
        print("Error happens when doing casename %s. Going to ignore this one..." % casenames[i])

        continue
    
    
    print("Dimension:", f.variables[args.varname_var])
    new_casenames.append([casenames[i], legends[i], colors[i], linestyles[i]])
    mean_datas.append((f.variables[args.varname_mean][indices] - args.y_offset) / args.yscale)
    std_datas.append(np.sqrt(f.variables[args.varname_var][indices]) / args.yscale)
    
    f.close()

casenames = new_casenames

f = Dataset(args.domain_file, "r")
lat = f.variables["yc"][:, 1]
f.close()

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12), sharex=True)

ax1.set_title("Mean of %s%s" % (args.display_varname, args.extra_title))
ax2.set_xlabel("Latitude [deg]")
ax2.set_xticks([-90, -60, -30, 0, 30, 60, 90])

ax1.set_ylabel(args.ylabel)
ax2.set_ylabel(args.ylabel)

for i, (casename, legend, color, linestyle) in enumerate(casenames): 
    ax1.plot(lat, mean_datas[i], linewidth=2, label=legend, color=color, linestyle=linestyle)
    ax2.plot(lat, std_datas[i],  linewidth=2, label=legend, color=color, linestyle=linestyle)

ax1.legend()

ax1.grid(True)
ax2.grid(True)

if args.y_range_mean != "":
    ax1.set_ylim(args.y_range_mean)

if args.y_range_std != "":
    ax2.set_ylim(args.y_range_std)


fig.savefig("%s/mc_meridional_mean_std_%s_%s%s.png" % (args.output_dir, args.varname_mean, args.varname_var, args.extra_filename), dpi=200)

#plt.show()
