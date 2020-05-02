import numpy as np
import matplotlib.pyplot as plt
import os

import matplotlib.style as style
style.use('ggplot')

from matplotlib import rc
rc({'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def read_out(filepath):
    seq_runtimes = list()
    top_down_runtimes = list()
    with open(filepath, "r") as f:
        lines = f.readlines()
        for line in lines:
            line_split = line.split(" ")
            data = [int(float(x)) for x in line_split]
            seq_runtimes.append(data[0])
            top_down_runtimes.append(data[1])
            print(data)
    seq_runtimes = np.array(seq_runtimes)
    top_down_runtimes = np.array(top_down_runtimes)
    # print("Average Top Down Runtime:", np.mean(top_down_runtimes)/1e6, "ms")
    # print("Average Bottom Up Runtime:", np.mean(bottom_up_runtimes)/1e6, "ms")
    # print("Average Hybrid Runtime:", np.mean(hybrid_runtimes)/1e6, "ms")
    return np.mean(seq_runtimes), np.mean(top_down_runtimes)

def plot_runtimes(prefix):
    ts = np.arange(1,9)
    avg_seq_runtimes = [590438225,590330997,590942967,590062792,589954669,589243878,591409536,592531934]
    avg_par_runtimes = [591566717,309022344,219178613,176232405,156674935,144091538,134204101,131053753]
    avg_seq_runtimes = list(map(lambda x : x/1e6, avg_seq_runtimes))
    avg_par_runtimes = list(map(lambda x : x/1e6, avg_par_runtimes))
    # for t in ts:
    #     avg_seq_runtime, avg_top_down_par_runtime = read_out(os.path.join(prefix, "le_lists_powerlaw2_" + str(t) + ".txt"))
    #     avg_seq_runtimes.append(avg_seq_runtime)
    #     avg_top_down_par_runtimes.append(avg_top_down_par_runtime)
    # plt.plot(ts, avg_seq_runtimes)
    plt.plot(ts, avg_seq_runtimes)
    plt.plot(ts, avg_par_runtimes)
    plt.xlabel("Number of Threads")
    plt.ylabel("Time (ms)")
    plt.legend(["Sequential", "Parallel"])
    plt.title("LE-Lists Runtimes on Random Powerlaw Graphs ($|V| = 20000$, $\mu_{deg} = 10$, $p = 0.25$)", fontsize=12)
    plt.savefig("le_lists_plots/le_lists_runtimes.png", dpi=500)
    plt.show()

prefix = "../results/le_lists/"
plot_runtimes(prefix)