import numpy as np
import matplotlib.pyplot as plt
import os

import matplotlib.style as style
style.use('ggplot')

from matplotlib import rc
rc({'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def read_out(filepath):
    top_down_runtimes = list()
    bottom_up_runtimes = list()
    hybrid_runtimes = list()

    top_down_seq_runtimes = list()
    top_down_par_runtimes = list()
    bottom_up_seq_runtimes = list()
    bottom_up_par_runtimes = list()
    hybrid_runtimes = list()

    with open(filepath, "r") as f:
        lines = f.readlines()
        for line in lines:
            line_split = line.split(" ")
            data = [int(float(x)) for x in line_split]
            top_down_seq_runtimes.append(data[0])
            top_down_par_runtimes.append(data[1])
            bottom_up_seq_runtimes.append(data[2])
            bottom_up_par_runtimes.append(data[3])
            hybrid_runtimes.append(data[4])
            # print(data)
    top_down_seq_runtimes = np.array(top_down_seq_runtimes)
    top_down_par_runtimes = np.array(top_down_par_runtimes)
    bottom_up_seq_runtimes = np.array(bottom_up_seq_runtimes)
    bottom_up_par_runtimes = np.array(bottom_up_par_runtimes)
    hybrid_runtimes = np.array(hybrid_runtimes)
    # print("Average Top Down Seq Runtime:", np.mean(top_down_seq_runtimes)/1e3, "microseconds")
    # print("Average Bottom Up Seq Runtime:", np.mean(bottom_up_seq_runtimes)/1e3, "microseconds")
    # print("Average Top Down Par Runtime:", np.mean(top_down_par_runtimes)/1e3, "microseconds")
    # print("Average Bottom Up Par Runtime:", np.mean(bottom_up_par_runtimes)/1e3, "microseconds")
    # print("Average Hybrid Runtime:", np.mean(hybrid_runtimes)/1e3, "microseconds")
    return np.mean(top_down_seq_runtimes), np.mean(top_down_par_runtimes), np.mean(bottom_up_seq_runtimes), np.mean(bottom_up_par_runtimes), np.mean(hybrid_runtimes)

def plot_speedup_internet(prefix):
    ts = np.arange(1,9)
    avg_top_down_seq_runtimes = list()
    avg_top_down_par_runtimes = list()
    avg_bottom_up_seq_runtimes = list()
    avg_bottom_up_par_runtimes = list()
    avg_hybrid_par_runtimes = list()
    for t in ts:
        avg_top_down_seq_runtime, avg_top_down_par_runtime, avg_bottom_up_seq_runtime, avg_bottom_up_par_runtime, avg_hybrid_par_runtime = read_out(os.path.join(prefix, "bfs_internet_" + str(t) + ".txt"))
        avg_top_down_seq_runtimes.append(avg_top_down_seq_runtime)
        avg_top_down_par_runtimes.append(avg_top_down_par_runtime)
        avg_bottom_up_seq_runtimes.append(avg_bottom_up_seq_runtime)
        avg_bottom_up_par_runtimes.append(avg_bottom_up_par_runtime)
        avg_hybrid_par_runtimes.append(avg_hybrid_par_runtime)
    plt.plot(ts, avg_top_down_seq_runtimes)
    plt.plot(ts, avg_top_down_par_runtimes)
    plt.plot(ts, avg_bottom_up_seq_runtimes)
    plt.plot(ts, avg_bottom_up_par_runtimes)
    plt.plot(ts, avg_hybrid_par_runtimes)
    plt.xlabel("Number of Threads")
    plt.ylabel("Time (ns)")
    plt.legend(["Top Down (Seq)", "Top Down (Par)", "Bottom Up (Seq)", "Bottom Up (Par)", "Hybrid (Par)"])
    plt.title("BFS Runtimes on Random Internet AS Graphs ($|V| = 10000$)")
    plt.savefig("bfs_plots/bfs_internet_runtimes.png", dpi=500)
    plt.show()

def plot_speedup_powerlaw(prefix):
    ts = np.arange(1,9)
    avg_top_down_seq_runtimes = list()
    avg_top_down_par_runtimes = list()
    avg_bottom_up_seq_runtimes = list()
    avg_bottom_up_par_runtimes = list()
    avg_hybrid_par_runtimes = list()
    for t in ts:
        avg_top_down_seq_runtime, avg_top_down_par_runtime, avg_bottom_up_seq_runtime, avg_bottom_up_par_runtime, avg_hybrid_par_runtime = read_out(os.path.join(prefix, "bfs_powerlaw_" + str(t) + ".txt"))
        avg_top_down_seq_runtimes.append(avg_top_down_seq_runtime)
        avg_top_down_par_runtimes.append(avg_top_down_par_runtime)
        avg_bottom_up_seq_runtimes.append(avg_bottom_up_seq_runtime)
        avg_bottom_up_par_runtimes.append(avg_bottom_up_par_runtime)
        avg_hybrid_par_runtimes.append(avg_hybrid_par_runtime)
    plt.plot(ts, avg_top_down_seq_runtimes)
    plt.plot(ts, avg_top_down_par_runtimes)
    plt.plot(ts, avg_bottom_up_seq_runtimes)
    plt.plot(ts, avg_bottom_up_par_runtimes)
    plt.plot(ts, avg_hybrid_par_runtimes)
    plt.xlabel("Number of Threads")
    plt.ylabel("Time (ns)")
    plt.legend(["Top Down (Seq)", "Top Down (Par)", "Bottom Up (Seq)", "Bottom Up (Par)", "Hybrid (Par)"])
    plt.title("BFS Runtimes on Random Powerlaw Graphs ($|V| = 20000$, $m = 5$, $p = 0.25$)", fontsize=12)
    plt.savefig("bfs_plots/bfs_powerlaw_runtimes.png", dpi=500)
    plt.show()

plot_speedup_internet("../results/bfs/")
# plot_speedup_powerlaw("../results/bfs/")
# read_out("../results/bfs/bfs.txt")