import os
import numpy as np
import matplotlib.pyplot as plt

results_dir = './timings'
time_results = {}
flops_results = {}
dev_results = {}
runs = 0

for timefile in os.listdir(results_dir):
    flops = []
    times = []
    testname, _= os.path.splitext(timefile)
    with open("%s/%s" % (results_dir, timefile), 'r') as f:
        f.readline()
        for line in f:
            time, flop, _ = line.split()
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
flops_bars = plt.bar(index + bar_width, flops_results.values(), bar_width,
                     alpha=opacity,
                     color='r', label='FLOPS')

dev_bars = plt.bar(index + bar_width, dev_results.values(), 0.001,
                   alpha=1,
                   color='b')


def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 1.01*height,
                '%.2f' % float(height),
                ha='center', va='bottom')


autolabel(time_bars)
autolabel(flops_bars)
autolabel(dev_bars)
plt.xticks(index + bar_width/2, time_results.keys(), rotation=80)
plt.xlabel('Transformation type')
plt.title('600 Grid, 500ms, Average performance over %d runs. gcc' % runs)
plt.legend()
plt.show()

