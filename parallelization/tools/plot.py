import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

if __name__ == "__main__":
    # data = pd.read_csv("out/foo.csv", sep=",")
    # data = np.array(data).transpose()
    # plt.plot(data[0], data[1], label="f(x)")

    # data = pd.read_csv("out/foo_dd_a.csv")
    # data = np.array(data).transpose()
    # plt.plot(data[0], data[1], label="f''(x) a")

    # data = pd.read_csv("out/foo_dd_n.csv")
    # data = np.array(data).transpose()
    # plt.plot(data[0], data[1], label="f''(x) n")

    plt.legend(loc="best")
    plt.show()

    for i in range(5, 12):
        data = pd.read_csv("out/num_conv_" + str(i) + ".csv")
        data = np.array(data).transpose()
        plt.plot(data[0], data[1], label=str(i))

    plt.legend(loc="best")
    plt.show()