import numpy as np
import matplotlib.pyplot as plt

files = [
    # 'MASJ0014_0_2016-08-22_r.calib',
    # 'MASJ0014_0_2016-08-23_r.calib',
    # 'MASJ0014_0_2016-08-24_r.calib',
    # 'MASJ0014_0_2016-08-25_r.calib',
    # 'MASJ0014_1_2016-11-07_r.calib',
    # 'MASJ0014_1_2017-06-09_r.calib',
    # 'MASJ0014_1_2017-06-11_r.calib'
]
binsize = len(files)


master_ts = []
master_fl = []
master_fe = []

fig, ax = plt.subplots()
for file in files:
    data = np.loadtxt(file, delimiter=' ')
    ts = data[:,0]
    fl = data[:,1]
    fe = data[:,2]

    ax.step(ts, fl, label=file)

    master_ts += list(ts)
    master_fl += list(fl)
    master_fe += list(fe)

plt.legend()
plt.show(block=False)

master_ts = np.array(master_ts)
master_fl = np.array(master_fl)
master_fe = np.array(master_fe)

sorted_args = np.argsort(master_ts)

master_ts = master_ts[sorted_args]
master_fl = master_fl[sorted_args]
master_fe = master_fe[sorted_args]


bin_ts = []
bin_fl = []
bin_fe = []

beg = 0
end = binsize
i = 0
while end < len(master_ts):
    beg = i*binsize
    i += 1
    end = i*binsize

    bin_ts.append(np.mean(master_ts[beg:end]))
    bin_fl.append(np.mean(master_fl[beg:end]))
    bin_fe.append(np.mean(master_fe[beg:end]))


fig, ax = plt.subplots()
ax.step(bin_ts, bin_fl)
plt.show(block=False)

cont = input("Write to a file? y/n: ")
if cont.lower()[0] == 'y':
    oname = input("Enter a filename: ")
    oname += '.calib'
    with open(oname, 'w') as f:
        f.write("# This file was produced by binning the following files by {}:\n".format(binsize))
        for cf in files:
            f.write("# {}\n".format(cf))
        f.write("#\n# phase, flux, error\n")
        for t, fl, fe in zip(bin_ts, bin_fl, bin_fe):
            f.write("{} {} {}\n".format(t, fl, fe))

plt.close()