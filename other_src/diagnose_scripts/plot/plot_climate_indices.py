import matplotlib.pyplot as plt
import Nio, sys, argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--data-file-PDO', dest='data_file_PDO')
parser.add_argument('--data-file-AO', dest='data_file_AO')
parser.add_argument('--output-dir', dest='output_dir')

args = parser.parse_args()

print(args)

f = Nio.open_file(args.data_file_PDO, "r")
PDO = f.variables["PDO"][:]
f.close()

f = Nio.open_file(args.data_file_AO, "r")
AO = f.variables["AO"][:]
f.close()


PDO /= np.std(PDO)
AO  /= np.std(AO)

fig, ax = plt.subplots(1, 1, figsize=(12, 8))


time = np.arange(len(PDO)) / 12
years = np.floor(len(time) / 12.0)

ax.set_title("PDO and AO indices (%d years)" % (years, ))

ax.plot([time[0], time[-1]], [0, 0], "-", color="#cccccc", linewidth=2)
ax.plot(time, PDO, "b-", linewidth=2, label="PDO")
ax.plot(time, AO, "r--", linewidth=2, label="AO")
ax.legend()

fig.savefig("%s/cimate_indices.png" % (args.output_dir,), dpi=200)
