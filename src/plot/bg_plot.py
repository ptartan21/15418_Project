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
    with open(filepath, "r") as f:
        lines = f.readlines()
        for line in lines:
            line_split = line.split(" ")
            data = [int(float(x)) for x in line_split]
            top_down_runtimes.append(data[0])
            bottom_up_runtimes.append(data[1])
            hybrid_runtimes.append(data[2])
            print(data)
    top_down_runtimes = np.array(top_down_runtimes)
    bottom_up_runtimes = np.array(bottom_up_runtimes)
    hybrid_runtimes = np.array(hybrid_runtimes)
    # print("Average Top Down Runtime:", np.mean(top_down_runtimes)/1e6, "ms")
    # print("Average Bottom Up Runtime:", np.mean(bottom_up_runtimes)/1e6, "ms")
    # print("Average Hybrid Runtime:", np.mean(hybrid_runtimes)/1e6, "ms")
    return np.mean(top_down_runtimes), np.mean(bottom_up_runtimes), np.mean(hybrid_runtimes)

def plot_runtimes_powerlaw2(prefix):
    ts = np.arange(1,9)
    # avg_seq_runtimes = list()
    # avg_top_down_par_runtimes = list()
    # avg_bottom_up_par_runtimes = list()
    # avg_hybrid_par_runtimes = list()
    # for t in ts:
    #     avg_top_down_par_runtime, avg_bottom_up_par_runtime, avg_hybrid_par_runtime = read_out(os.path.join(prefix, "bg_powerlaw2_" + str(t) + ".txt"))
    #     # avg_seq_runtimes.append(avg_seq_runtime)
    #     avg_top_down_par_runtimes.append(avg_top_down_par_runtime)
    #     avg_bottom_up_par_runtimes.append(avg_bottom_up_par_runtime)
    #     avg_hybrid_par_runtimes.append(avg_hybrid_par_runtime)
    avg_bottom_up_par_runtimes = [8448404.000000,5384894.000000,4335959.000000,4044581.000000,3553171.000000,3327946.000000,3193336.000000,2868608.000000]
    avg_bottom_up_par_runtimes = list(map(lambda x : x/1e6, avg_bottom_up_par_runtimes))
    # plt.plot(ts, avg_seq_runtimes)
    # plt.plot(ts, avg_top_down_par_runtimes)
    plt.plot(ts, avg_bottom_up_par_runtimes)
    # plt.plot(ts, avg_hybrid_par_runtimes)
    plt.xlabel("Number of Threads")
    plt.ylabel("Time (ms)")
    plt.legend(["Bottom Up (Par)"])
    plt.title("Ball Growing Runtimes on Random Powerlaw Graphs ($|V| = 20000$, $\mu_{deg} = 10$, $p = 0.25$)", fontsize=12)
    plt.savefig("bg_plots/bg_powerlaw2_runtimes.png", dpi=500)
    plt.show()

def plot_runtimes_random(prefix):
    ts = np.arange(1,9)
    avg_seq_runtimes = list()
    avg_top_down_par_runtimes = list()
    avg_bottom_up_par_runtimes = list()
    avg_hybrid_par_runtimes = list()
    for t in ts:
        avg_top_down_par_runtime, avg_bottom_up_par_runtime, avg_hybrid_par_runtime = read_out(os.path.join(prefix, "bg_random_" + str(t) + ".txt"))
        # avg_seq_runtimes.append(avg_seq_runtime)
        avg_top_down_par_runtimes.append(avg_top_down_par_runtime)
        avg_bottom_up_par_runtimes.append(avg_bottom_up_par_runtime)
        avg_hybrid_par_runtimes.append(avg_hybrid_par_runtime)
    # plt.plot(ts, avg_seq_runtimes)
    plt.plot(ts, avg_top_down_par_runtimes)
    plt.plot(ts, avg_bottom_up_par_runtimes)
    plt.plot(ts, avg_hybrid_par_runtimes)
    plt.xlabel("Number of Threads")
    plt.ylabel("Time (ns)")
    plt.legend(["Top Down (Par)", "Bottom Up (Par)", "Hybrid (Par)"])
    plt.title("Ball Growing Runtimes on Random Graph ($|V| = 20000$, $|E| = 100000$)")
    plt.savefig("bg_plots/bg_random_runtimes.png", dpi=500)
    plt.show()

prefix = "../results/bg/"
plot_runtimes_powerlaw2(prefix)
# plot_runtimes_random(prefix)