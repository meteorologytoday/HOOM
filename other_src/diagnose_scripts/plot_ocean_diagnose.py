import Ngl, Nio
import sys, argparse
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument('--data-file-SSTAYYC', dest='SSTAYYC_file')
parser.add_argument('--data-file-SSTAVAR', dest='SSTAVAR_file')
parser.add_argument('--domain-file', dest='domain_file')
parser.add_argument('--output-dir', dest='output_dir')

args = parser.parse_args()

print(str(args))


g = Nio.open_file(args.domain_file, "r")

lon = g.variables["xc"][1, :]                   #-- read clon
lat = g.variables["yc"][:, 1]                   #-- read clat


f = Nio.open_file(args.SSTAYYC_file, "r")
data = f.variables["SSTAYYC"][1, :, :]
missing_value = f.variables["SSTAYYC"]._FillValue[0]
data[np.isnan(data)] = missing_value
f.close()

print("Creating workstation...")

wks_type = "png"
wks      =  Ngl.open_wks(wks_type, "%s/SSTAYYC" % (args.output_dir,)) #-- open a workstation

print("Defining res...")
cnres                 = Ngl.Resources()

cnres.sfMissingValueV = missing_value

# Contour resources
cnres.cnFillOn        = True
cnres.cnFillPalette   = "BlueYellowRed"      # New in PyNGL 1.5.0
cnres.cnLinesOn       = False
cnres.cnLineLabelsOn  = False

# Labelbar resource
cnres.lbOrientation   = "horizontal"

# Scalar field resources
cnres.sfXArray        = lon
cnres.sfYArray        = lat

# Map resources
cnres.mpFillOn               = True
cnres.mpFillDrawOrder        = "PostDraw"
cnres.mpLandFillColor        = "Gray"
cnres.mpOceanFillColor       = "Transparent"
cnres.mpInlandWaterFillColor = "Transparent"

print("Start plotting...")
contour = Ngl.contour_map(wks, data, cnres)
#Ngl.end()
print("plotting done.")



# ================================================


f = Nio.open_file(args.SSTAVAR_file, "r")
data = f.variables["SSTAVAR"][1, :, :]
missing_value = f.variables["SSTAVAR"]._FillValue[0]
data[np.isnan(data)] = missing_value
f.close()

print("Creating workstation...")

wks_type = "png"
wks      =  Ngl.open_wks(wks_type, "%s/SSTAVAR" % (args.output_dir,))

print("Defining res...")
cnres                 = Ngl.Resources()

cnres.sfMissingValueV = missing_value

# Contour resources
cnres.cnFillOn        = True
cnres.cnFillPalette   = "BlueYellowRed"      # New in PyNGL 1.5.0
cnres.cnLinesOn       = False
cnres.cnLineLabelsOn  = False

# Labelbar resource
cnres.lbOrientation   = "horizontal"

# Scalar field resources
cnres.sfXArray        = lon
cnres.sfYArray        = lat

# Map resources
cnres.mpFillOn               = True
cnres.mpFillDrawOrder        = "PostDraw"
cnres.mpLandFillColor        = "Gray"
cnres.mpOceanFillColor       = "Transparent"
cnres.mpInlandWaterFillColor = "Transparent"

print("Start plotting...")
contour = Ngl.contour_map(wks, data, cnres)
Ngl.end()
print("plotting done.")
