import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import collections

resultsdir = sys.argv[1]

def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.01*height,
                '%.3f' % float(height),
                ha='center', va='bottom')

# for timedir in os.listdir(resultsdir):
time_results = {}
flops_results = {}
dev_results = {}
runs = 0
timedir=""
for timefile in os.listdir("%s/%s" %(resultsdir, timedir)):
    flops = []
    times = []
    testname, _= os.path.splitext(timefile)
    with open("%s/%s/%s" % (resultsdir, timedir, timefile), 'r') as f:
        for line in f:
            time, flop = line.split()
            times.append(float(time))
            flops.append(float(flop))
    if len(times) == 0:
        continue
    runs = len(times)
    np_time = np.array(times)
    meantime = np.mean(np_time)
    flops_results[meantime] = np.mean(np.array(flops))
    time_results[meantime] = testname
    dev_results[meantime] = np.std(np_time)

flops_results = collections.OrderedDict(sorted(flops_results.items(), reverse=True))
time_results = collections.OrderedDict(sorted(time_results.items(), reverse=True))
dev_results = collections.OrderedDict(sorted(dev_results.items(), reverse=True))


index = np.arange(len(time_results))
fig, ax = plt.subplots()
bar_width = 0.35
opacity = 0.5

keys = sorted(time_results, key=time_results.get, reverse=True)

time_bars = plt.bar(index, time_results.keys(), bar_width,
                    alpha=opacity,
                    color='g', label='Runtime (s)')
dev_bars = plt.bar(index + bar_width, dev_results.values(), 0.001,
                    alpha=1,
                    color='b')
flops_bars = plt.bar(index + bar_width, flops_results.values(), bar_width,
                    alpha=opacity,
                    color='r', label='GFLOPS')


autolabel(time_bars)
autolabel(flops_bars)
autolabel(dev_bars)
plt.xticks(index + bar_width/2, time_results.values(), rotation=80)
plt.xlabel('Transformation type')
plt.title('Average performance over %d runs. Compiled with icc' % runs)
plt.legend()
plt.subplots_adjust(bottom=0.27, top=0.95, right=0.95, left=0.05)
plt.show()
# plt.savefig('results.png')
print("##STDEV###")
print (dev_results)

