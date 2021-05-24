import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def plot_Q6_standard():
    data = np.genfromtxt("out/RK.csv", delimiter=",", skip_header=1, comments="#").transpose()

    start = data[0][0]
    stop = data[0][-1]
    x = np.arange(start, stop, 0.1)
    y_a = np.cos(x)
    plt.plot(data[0], data[1], label="num")
    plt.plot(x, y_a, label="ana", linestyle="dashed")

    plt.xlabel("t")
    plt.ylabel("q(t)")

    plt.title("dt=0.5, m=1")

    plt.legend(loc="best")
    plt.show()


if __name__ == "__main__":
    plot_Q6_standard()

    # for i in range(5, 12):
    #     data = pd.read_csv("out/num_conv_" + str(i) + ".csv")
    #     data = np.array(data).transpose()
    #     plt.plot(data[0], data[1], label=str(i))

    # plt.legend(loc="best")
    # plt.show()