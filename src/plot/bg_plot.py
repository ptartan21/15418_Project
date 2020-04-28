import numpy as np
import matplotlib.pyplot as plt

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
    print("Average Top Down Runtime:", np.mean(top_down_runtimes)/1e6, "ms")
    print("Average Bottom Up Runtime:", np.mean(bottom_up_runtimes)/1e6, "ms")
    print("Average Hybrid Runtime:", np.mean(hybrid_runtimes)/1e6, "ms")

read_out("../results/ball_growing/bg.txt")