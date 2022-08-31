"""
Plot the Head and Tail.

Input the file 'n_b' consisting of the mean and std of r_1,...r_{n-1} of n-dim BKZ_b bases
"""

import matplotlib.pyplot as plt
from math import sqrt


def stat(l):
    """Return mean and standard deviation of elements in list l."""
    ans = 0.0
    ave = sum(l) / len(l)
    for x in l:
        ans += (x - ave)**2
    ans /= len(l)
    return (ave, sqrt(ans))


def marker(i):
    """Different markers."""
    if i == 0:
        return "^"
    if i == 1:
        return "s"
    if i == 2:
        return "v"
    if i == 3:
        return "o"
    if i == 4:
        return ">"
    if i == 5:
        return "D"
    if i == 6:
        return "<"
    if i == 7:
        return "p"


def plot_ht(l_n, l_beta ,gen ):
    """Plot head and tail."""
    MAX = max(l_n)
    for beta in l_beta:
        fig = plt.figure(figsize=(8, 5))
        if gen == "gen_ntru_instance_circulant":
            plt.title("Head and Tail in " + "BKZ_" + str(beta) + "ntru_circulant" +  " basis", fontsize=20)
        elif gen == "gen_ntru_instance_matrix":
            plt.title("Head and Tail in " + "BKZ_" + str(beta) + "ntru_matrix" +  " basis", fontsize=20)
        elif gen == "uniform":
            plt.title("Head and Tail in " + "BKZ_" + str(beta) + "ntru_uniform " +  " basis", fontsize=20)
        i = 0
        for n in l_n:
            h_Av = []
            h_Std = []
            t_Av = []
            t_Std = []
            if gen == "gen_ntru_instance_circulant":
                f = open("stat_r/" + str(n) +"/ntru_circulant/"+str(n) + "_" + str(beta) )
            elif gen == "gen_ntru_instance_matrix":
                f = open("stat_r/" + str(n) +"/ntru_matrix/"+str(n) + "_" + str(beta) )
            elif gen == "uniform":
                f = open("stat_r/" + str(n) +"/ntru_uniform/"+str(n) + "_" + str(beta) )
            cnt = 0
            for line in f:
                if line == "":
                    break
                data = line.split()
                print(data[0])
                if cnt < n / 2 - 1:
                    h_Av += [float(data[0])]
                    h_Std += [float(data[1])]
                elif cnt >= n / 2 - 1:
                    t_Av += [float(data[0])]
                    t_Std += [float(data[1])]
                cnt += 1
            plt.plot(range(1, int (n / 2) ), h_Std, "b" + marker(i + 4) + "-", markersize=5, linewidth=1)
            plt.plot(range(MAX - 1 - n + int(n / 2) + 1, MAX), t_Std, "b" + marker(i + 4) + ":", markersize=5, linewidth=1, label='Std for n=' + str(n))
            i += 1
        #plt.legend(loc=1, ncol=2, prop={'size': 16})
        #plt.ylim(-0.02, 0.14)
        #plt.plot([beta, beta], [0, 0.125], "k--")
        #plt.plot([MAX - beta, MAX - beta], [0, 0.125], "k--")
        #plt.xticks([0, 20, 40, 60, 80, 100, 120, 140])
        #ax = fig.add_subplot(1, 1, 1)
        #xlabels = [0, 20, 40,  n-40, n-20, n]
        #label_format = '{:,.1f}'
        #ax.set_xticklabels([label_format.format(x) for x in xlabels])
        #ax.set_xticklabels([0, 20, 40, 60, "n-60", "n-40", "n-20", "n"])
        #ax.set_xticklabels([0, 20, 40,  n-40, n-20, n])
        if gen == "gen_ntru_instance_circulant":
            plt.savefig("fig_ht/ht_in_bkz_" + str(beta) +"ntru_circulant "  + '.png')
        elif gen == "gen_ntru_instance_matrix":
            plt.savefig("fig_ht/ht_in_bkz_" + str(beta) +"ntru_matrix "  + '.png')
        elif gen == "uniform":
            plt.savefig("fig_ht/ht_in_bkz_" + str(beta) +"ntru_uniform "  + '.png')
        plt.close()


def my_plot_ht(l_n, l_beta ,gen ):
    """Plot head and tail."""
    MAX = max(l_n)
    for beta in l_beta:
        fig = plt.figure(figsize=(8, 5))
        if gen == "gen_ntru_instance_circulant":
            plt.title("Head and Tail in " + "BKZ_" + str(beta) + "ntru_circulant" +  " basis", fontsize=20)
        elif gen == "gen_ntru_instance_matrix":
            plt.title("Head and Tail in " + "BKZ_" + str(beta) + "ntru_matrix" +  " basis", fontsize=20)
        elif gen == "uniform":
            plt.title("Head and Tail in " + "BKZ_" + str(beta) + "ntru_uniform " +  " basis", fontsize=20)
        i = 0
        for n in l_n:
            h_Av = []
            h_Std = []
            t_Av = []
            t_Std = []
            if gen == "gen_ntru_instance_circulant":
                f = open("stat_r/" + str(n) +"/ntru_circulant/"+str(n) + "_" + str(beta) )
            elif gen == "gen_ntru_instance_matrix":
                f = open("stat_r/" + str(n) +"/ntru_matrix/"+str(n) + "_" + str(beta) )
            elif gen == "uniform":
                f = open("stat_r/" + str(n) +"/ntru_uniform/"+str(n) + "_" + str(beta) )
            cnt = 0
            for line in f:
                if line == "":
                    break
                data = line.split()
                #print(data[0])
                if cnt < n / 2 - 1:
                    h_Av += [float(data[0])]
                    h_Std += [float(data[1])]
                elif cnt >= n / 2 - 1:
                    t_Av += [float(data[0])]
                    t_Std += [float(data[1])]
                cnt += 1
            plt.ylim(-0.05, 0.16)
            plt.plot(range(1, int (n / 2) ), h_Av,     '>:r', markersize=5, linewidth=1)
            plt.plot(range(MAX - 1 - n + int(n / 2) + 1, MAX), t_Av,  '>:r', markersize=5, linewidth=1, label='Std for n=' + str(n))
            plt.plot(range(1, int (n / 2) ), h_Std, '<:b', markersize=5, linewidth=1)
            plt.plot(range(MAX - 1 - n + int(n / 2) + 1, MAX), t_Std, '<:b', markersize=5, linewidth=1, label='Std for n=' + str(n))
            
            i += 1
        if gen == "gen_ntru_instance_circulant":
            plt.savefig("fig_ht/140ht_in_bkz_" + str(beta) +"ntru_circulant "  + '.png')
        elif gen == "gen_ntru_instance_matrix":
            plt.savefig("fig_ht/140ht_in_bkz_" + str(beta) +"ntru_matrix "  + '.png')
        elif gen == "uniform":
            plt.savefig("fig_ht/140ht_in_bkz_" + str(beta) +"ntru_uniform "  + '.png')
        plt.close()


N = [140]
#N = [140,  100]
Beta = [5,10,15,20,25,30]
#my_plot_ht(N, Beta , gen = "gen_ntru_instance_circulant" )
#my_plot_ht(N, Beta , gen = "gen_ntru_instance_matrix" )
my_plot_ht(N, Beta , gen = "uniform" )