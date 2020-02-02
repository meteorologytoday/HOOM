import cartopy.crs as ccrs
import matplotlib as mplt
mplt.use('Agg')

import matplotlib.pyplot as plt
from matplotlib import cm

from netCDF4 import Dataset

import sys, argparse
import numpy as np

def ext(data):
    s = data.shape
    ndata = np.zeros((s[0], s[1]+1))
    ndata[:, 0:-1] = data
    ndata[:, -1] = data[:, 0]
    return ndata
 

def ext_axis(lon):
    return np.append(lon, 360) 
 
parser = argparse.ArgumentParser()

parser.add_argument('--data-file')
parser.add_argument('--domain-file')
parser.add_argument('--output-dir')
parser.add_argument('--output-file-prefix', type=str, default="")
parser.add_argument('--title-prefix', type=str, default="")

args = parser.parse_args()

print(args)

if args.output_file_prefix != "":
    args.output_file_prefix = "%s_" % args.output_file_prefix


with Dataset(args.domain_file, "r") as g:
    lon = g.variables["xc"][1, :]                   #-- read clon
    lat = g.variables["yc"][:, 1]                   #-- read clat

    lon = ext_axis(lon)

with Dataset(args.data_file, "r") as f:
    data = f.variables["ice_cov"][:] * 100.0

# I think this is a bug in cartopy that projection are not consistent
proj1 = ccrs.PlateCarree(central_longitude=210.0)
proj2 = ccrs.PlateCarree(central_longitude=0.0)

cmap = cm.get_cmap("YlGnBu")
cmap.set_under("white")


clevels = np.linspace(15, 100, 18)


months = "Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec".split(" ")

for m in range(len(months)):
   
    print("Doing month: %d" % (m+1,)) 
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6), subplot_kw={'projection': proj1, 'aspect': 1.5})
    ax.set_title("%sSea-ice of month %s" % (args.title_prefix, months[m]))

#    print(data.shape)
 #   print(lon.shape)
  #  print(lat.shape)

    mappable = ax.contourf(lon, lat, ext(data[m, :, :]), clevels, cmap=cmap, extend="both", transform=proj2)

    ax.set_xlim([-90, 90])
    ax.set_xticks(np.linspace(-180, 180, 13))
    ax.set_xticklabels(["%d" % ((d+30)%360,) for d in np.linspace(0, 360, 13)])

    ax.set_yticks(np.linspace(-90,  90, 7))
    ax.set_yticklabels(["%d" % d for d in np.linspace(-90,  90, 7)])

    ax.set_global()
    ax.coastlines()
    
    cb = fig.colorbar(mappable, ax=ax, orientation="vertical", ticks=clevels)
    cb.ax.set_ylabel('Sea-ice fraction [$ \%  $]', rotation=90)

    filename = "%s/%sSI-%02d.png" % (args.output_dir, args.output_file_prefix, m+1)
    fig.savefig(filename, dpi=300, transparent=True)
    plt.close(fig)


# Plot mean

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 6), subplot_kw={'projection': proj1, 'aspect': 1.5})
ax.set_title("%s Mean Sea-ice" % (args.title_prefix,))

mappable = ax.contourf(lon, lat, ext(data.mean(axis=0)), clevels, cmap=cmap, extend="both", transform=proj2)

ax.set_xlim([-90, 90])
ax.set_xticks(np.linspace(-180, 180, 13))
ax.set_xticklabels(["%d" % ((d+30)%360,) for d in np.linspace(0, 360, 13)])

ax.set_yticks(np.linspace(-90,  90, 7))
ax.set_yticklabels(["%d" % d for d in np.linspace(-90,  90, 7)])

ax.set_global()
ax.coastlines()

cb = fig.colorbar(mappable, ax=ax, orientation="vertical", ticks=clevels)
cb.ax.set_ylabel('Sea-ice fraction [$ \%  $]', rotation=90)

filename = "%s/%sSI-mean.png" % (args.output_dir, args.output_file_prefix)
fig.savefig(filename, dpi=300, transparent=True)
plt.close(fig)
