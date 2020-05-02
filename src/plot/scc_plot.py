# hey ;)
import matplotlib.pyplot as plt

import matplotlib.style as style
style.use('ggplot')

from matplotlib import rc
rc({'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def plot_scc_runtimes():
    num_threads = [1,2,3,4,5,6,7,8]
    seq_bu = [14071002, 14082249, 14146934, 14277245, 14390688, 14520803, 14497715, 14499155]
    seq_td = [7139864, 7127555, 7083515, 7083284, 7085122, 7072155, 7036239, 7068205]
    par_bu = [14273671, 8581526, 6711274, 5710824, 5147140, 4768797, 4552728, 4281837]
    par_td = [10305774, 7422807, 6351799, 5452563, 4916161, 4583268, 4413883, 4303516]
    hybrid = [8342077, 5643662, 4606685, 3343321, 3261011, 3104153, 3112135, 3060815]

    seq_bu = list(map(lambda x : x/1e6, seq_bu))
    seq_td = list(map(lambda x : x/1e6, seq_td))
    par_bu = list(map(lambda x : x/1e6, par_bu))
    par_td = list(map(lambda x : x/1e6, par_td))
    hybrid = list(map(lambda x : x/1e6, hybrid))


    plt.title('SCC-Detection Runtimes on Watts-Strogatz Graphs ($|V| = 20000$, $k = 20$, $p = 0.25$)', fontsize=12) 
    plt.ylabel('Time (ms)')
    plt.xlabel('Number of Threads')
    b, = plt.plot(num_threads, seq_bu)
    a, = plt.plot(num_threads, seq_td)
    d, = plt.plot(num_threads, par_bu)
    c, = plt.plot(num_threads, par_td)
    e, = plt.plot(num_threads, hybrid)
    plt.legend([b,a,d,c,e], ['Bottom Up (Seq)', 'Top Down (Seq)', 'Bottom Up (Par)', 'Top Down (Par)', 'Hybrid (Par)'], loc = 'upper right')
    plt.savefig("scc_plots/scc_ws_runtimes.png", dpi=500)
    plt.show()

def plot_scc_runtimes_small():
    num_threads = [1,2,3,4,5,6,7,8]
    seq_bu = [3022212, 3044468, 3068467, 3105594, 3196143, 3181999, 3174073, 3173763]
    seq_td = [1669195, 1702963, 1715536, 1725145, 1773203, 1765494, 1760036, 1763477] 
    par_bu = [3128788, 2235705, 1830850, 1592401, 1561855, 1475017, 1416089, 1358576] 
    par_td = [2466597, 2049554, 1716849, 1501517, 1398425, 1365566, 1304692, 1289207]
    hybrid = [1877913, 1730596, 1530764, 1337517, 1267433, 1184103, 1167146, 1157360]

    seq_bu = list(map(lambda x : x/1e6, seq_bu))
    seq_td = list(map(lambda x : x/1e6, seq_td))
    par_bu = list(map(lambda x : x/1e6, par_bu))
    par_td = list(map(lambda x : x/1e6, par_td))
    hybrid = list(map(lambda x : x/1e6, hybrid))

    plt.title('SCC-Detection Runtimes on Watts-Strogatz Graphs ($|V| = 5000$, $k = 20$, $p = 0.25$)', fontsize=12) 
    plt.ylabel('Time (ms)')
    plt.xlabel('Number of Threads')
    b, = plt.plot(num_threads, seq_bu)
    a, = plt.plot(num_threads, seq_td)
    d, = plt.plot(num_threads, par_bu)
    c, = plt.plot(num_threads, par_td)
    e, = plt.plot(num_threads, hybrid)
    plt.legend([b,a,d,c,e], ['Bottom Up (Seq)', 'Top Down (Seq)', 'Bottom Up (Par)', 'Top Down (Par)', 'Hybrid (Par)'], loc = 'upper right')
    plt.savefig("scc_plots/scc_ws_small_runtimes.png", dpi=500)
    plt.show()

plot_scc_runtimes()
plot_scc_runtimes_small()