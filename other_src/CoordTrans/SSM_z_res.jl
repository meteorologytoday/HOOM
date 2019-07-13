
# This z resolution is get from
zs_NCAR_LENS = - [
    0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 
    11000, 12000, 13000, 14000, 15000, 16000, 17019.681640625, 
    18076.12890625, 19182.125, 20349.931640625, 21592.34375, 22923.3125, 
    24358.453125, 25915.580078125, 27615.259765625, 29481.470703125, 
    31542.373046875, 33831.2265625, 36387.47265625, 39258.046875, 
    42498.88671875, 46176.65625, 50370.6875, 55174.91015625, 60699.66796875, 
    67072.859375, 74439.8046875, 82960.6953125, 92804.3515625, 
    104136.8203125, 117104.015625, 131809.359375, 148290.078125, 
    166499.203125, 186301.4375, 207487.390625, 229803.90625, 252990.40625, 
    276809.84375, 301067.0625, 325613.84375, 350344.875, 375189.1875, 
    400101.15625, 425052.46875, 450026.0625, 475012, 500004.6875, 
    525000.9375, 549999.0625
] / 100.0;


#############################################################################
# 2019/05/09
# I want to double the resultion below 2000 meters

# zs_NCAR_LENS[45] = 1863.0
# zs_NCAR_LENS[46] = 2074.9
cut = 45

zs_SSM = []
for i = 1:cut
    append!(zs_SSM, zs_NCAR_LENS[i])
    append!(zs_SSM, (zs_NCAR_LENS[i] + zs_NCAR_LENS[i+1]) / 2.0)
end

for i = cut+1:length(zs_NCAR_LENS)
    append!(zs_SSM, zs_NCAR_LENS[i])
end

#############################################################################

# 2019/05/10
# Make 40 layers (41 points). 5 meters increament in top 100m (20 layers),
# and a linearly increasing increament in the next 20 layers until 1000m.
zs_SSM = zeros(Float64, 41)
zs_SSM[1] = 0.0

dh = zeros(Float64, 40)
dh[ 1:20] .= 5.0
dh[21:40] = 5.0 .+ collect(1:20) * 80.0/21.0
for i = 1:length(dh)
    zs_SSM[i+1] = zs_SSM[i] - dh[i]
end


#############################################################################

# 2019/07/11
# Make 40 layers (41 points). Simply use the same zs as LENS (top 41 pts).

zs_SSM[:] = zs_NCAR_LENS[1:41] 

#############################################################################


zs_mid_NCAR_LENS = (zs_NCAR_LENS[2:end] + zs_NCAR_LENS[1:end-1] ) / 2.0
zs_mid_SSM = (zs_SSM[2:end] + zs_SSM[1:end-1] ) / 2.0




