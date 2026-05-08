import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import Counter

def read_file(path):
    data = {}

    with open(path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    current_name = None
    current_mode = None
    current_function = None
    table_start = False  # Indicates when the table data begins

    for line in lines:
        line = line.strip()
        if line.startswith('Name:'):
            # Update the current name
            path_name = line.split(':', 1)[1].strip()
            current_name = (os.path.dirname(path_name)).split('/')[-1]
            if current_name not in data:
                data[current_name] = {}
        elif line.startswith('Mode:'):
            # Update the current mode
            current_mode = line.split(':', 1)[1].strip()
            if current_name and current_mode not in data[current_name]:
                data[current_name][current_mode] = {}
        elif line.startswith('Function:'):
            # Update the current function
            current_function = line.split(':', 1)[1].strip()
            if current_name and current_mode and current_function not in data[current_name][current_mode]:
                data[current_name][current_mode][current_function] = {
                    "Ranges": [],
                    "Start": [],
                    "Finish": [],
                    "T2s": []
                }
        elif line == '':
            # The table data starts after this line
            table_start = True
        elif table_start and line:
            # Parse table data
            parts = line.split('\t')
            if len(parts) >= 3:
                start = float(parts[0])
                finish = float(parts[1])
                t2 = float(parts[2])
                data[current_name][current_mode][current_function]["Ranges"].append((start, finish))
                data[current_name][current_mode][current_function]["Start"].append(start)
                data[current_name][current_mode][current_function]["Finish"].append(finish)
                data[current_name][current_mode][current_function]["T2s"].append(t2)

    return data

def filter_values(T2, range, start, finish):
    filtered_T2     = []
    filtered_range  = []
    filtered_start  = []
    filtered_finish = []

    for i, tau in enumerate(T2):
        if 6 < tau < 50:
            filtered_T2.append(tau)
            filtered_range.append(range[i])
            filtered_start.append(start[i])
            filtered_finish.append(finish[i])

    return filtered_T2, filtered_range, filtered_start, filtered_finish

def time_analysis():
    ii = 0
    jj = 0
    fig_h, axes_h = plt.subplots(3, 3, figsize=(13, 7),sharex='col')
    fig_cm, axes_cm = plt.subplots(3, 3, figsize=(13, 7))

    for sample in data.keys():

        counts_ranges = []
        alpha = 1
        i = 0

        for mode in data[sample].keys():
            for function in data[sample][mode].keys():

                number_of_data = len(data[sample][mode].keys()) * 2
                colors = mpl.colormaps.get_cmap('jet').resampled(number_of_data)
                transpar = alpha - i / (number_of_data * 1.5)

                Ranges = data[sample][mode][function]["Ranges"]
                T2 = data[sample][mode][function]["T2s"]
                Start = data[sample][mode][function]["Start"]
                Finish = data[sample][mode][function]["Finish"]

                filtered_T2, filtered_range, filtered_start, filtered_finish = filter_values(T2, Ranges, Start, Finish)

                unique_values, ind, counts = np.unique(
                    np.round(filtered_T2, decimals=1), return_counts=True, return_index=True
                )

                hist_ax = axes_h[ii, jj]

                hist_ax.hist(
                    filtered_T2,
                    bins=len(counts),
                    color=colors(i),
                    label=f"{mode}: {function}",
                    alpha=transpar,
                    ec='k',
                    histtype='stepfilled'
                )

                i += 1

                for l, n in zip(ind, counts):
                    counts_ranges.extend([filtered_range[l]] * n)


        hist_ax.set_title(f'{sample}')


        # Heatmap
        unique_ranges = []
        counts = []
        unique_ranges, counts = np.unique(counts_ranges, axis=0, return_counts=True)

        start_vals = np.unique(unique_ranges[:, 0])
        finish_vals = np.unique(unique_ranges[:, 1])
        heatmap = np.zeros((len(start_vals), len(finish_vals)))

        for (start, finish), count in zip(unique_ranges, counts):
            start_idx = np.where(start_vals == start)[0][0]
            finish_idx = np.where(finish_vals == finish)[0][0]
            heatmap[start_idx, finish_idx] = count

        X, Y = np.meshgrid(finish_vals, start_vals)

        heat_axes = axes_cm[ii, jj]
        heatmap_plot = heat_axes.pcolormesh(X, Y, heatmap, cmap='jet', shading='auto')
        plt.colorbar(heatmap_plot, ax=heat_axes, label='Counts')
        heat_axes.set_title(f'{sample}')

        ii+=1
        if ii > 2:
            jj+=1
            ii = 0

    axes_h[1,2].legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig_h.supxlabel('T₂*, μs')
    fig_h.supylabel('Counts')
    fig_h.tight_layout()

    fig_cm.supxlabel('Finish')
    fig_cm.supylabel('Start')
    fig_cm.tight_layout()


# path = "statistical_analysis.txt"
# data = read_file(path)
# time_analysis()

# path2 = "statistical_analysis_full.txt"
# data = read_file(path2)
# time_analysis()

path1 = "statistical_analysis.txt"
data = read_file(path1)
time_analysis()
plt.show()

# add here total counts
print('done')