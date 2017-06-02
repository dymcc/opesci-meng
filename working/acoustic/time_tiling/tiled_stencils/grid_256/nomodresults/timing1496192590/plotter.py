import os
import numpy as np
import matplotlib.pyplot as plt

resultsdir = './timings'

def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.01*height,
                '%.3f' % float(height),
                ha='center', va='bottom')

time_results = {}
flops_results = {}
dev_results = {}
runs = 0
for timefile in os.listdir(resultsdir):
    flops = []
    times = []
    testname, _= os.path.splitext(timefile)
    with open("%s/%s" % (resultsdir, timefile), 'r') as f:
        for line in f:
            time, flop = line.split()
            times.append(float(time))
            flops.append(float(flop))
    if len(times) == 0:
        continue
    runs = len(times)
    np_flop = np.array(flops)
    np_time = np.array(times)
    flops_results[testname] = np.mean(np_flop)
    time_results[testname] = np.mean(np_time)
    dev_results[testname] = np.std(np_time)

print(flops_results)
print(time_results)

index = np.arange(len(time_results))
fig, ax = plt.subplots()
bar_width = 0.35
opacity = 0.5

time_bars = plt.bar(index, time_results.values(), bar_width,
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
plt.xticks(index + bar_width/2, time_results.keys(), rotation=80)
plt.xlabel('Transformation type')
plt.title('Average performance over %d runs. Compiled with icc' % runs)
plt.legend()
plt.show()
print("##STDEV###")
print (dev_results)

