import numpy as np
import matplotlib
#cols = 180 # number of column  =  the dimension of lattice
#divided_ch = '  ' # divided_character between numbers
from math import sqrt, log
import time
import matplotlib.pyplot as plt


def stat(l):
    """Return mean and standard deviation of elements in list l."""
    ans = 0.0
    ave = sum(l) / len(l)
    for x in l:
        ans += (x - ave)**2
    ans /= len(l)
    return (ave, sqrt(ans))




def dat_to_matrix(filename,cols , divided_ch = '  '):
    file = open(filename)
    lines = file.readlines()
    rows = len(lines)  
    datamat = np.zeros((rows, cols)) # 0:rows-1
    row = 0
    for line in lines:  
        line = line.strip().split(divided_ch) # strip remove block space in line
        datamat[row, :] = line[:]
        row += 1
        #print(len(line))  = n-1
    return datamat

def dat_to_one_list(filename, cols , divided_ch = '  '):
    file = open(filename)
    data_one = np.zeros(cols,dtype = float)
    lines = file.readlines()
    data_one = lines[0].strip().split(divided_ch)
    return data_one

def gen_cov(n, beta, i, j, q,gen,times ):
    """Return covariance between r_i and r_j."""
    if gen == "gen_ntru_instance_circulant":
        #f = open("stat_r/" + str(n) +"/ntru_circulant/"+str(n) + "_" + str(beta), "w")
        filename = "data-bkz1/"+"dimension" +str(n)+ "/" + "ntru_circulant/" + "r_i/"+str(times)+"times" +str(n)+ "n" +"_" + "q" + str(q) + "_"+ "beta"+ str(beta)  
    elif gen == "gen_ntru_instance_matrix":
        #f = open("stat_r/" + str(n) +"/ntru_matrix/"+str(n) + "_" + str(beta), "w")
        filename =  "data-bkz1/"+"dimension" +str(n)+ "/" + "ntru_matrix/" + "r_i/"+str(times)+"times" +str(n)+ "n" +"_" + "q" + str(q) + "_"+ "beta"+ str(beta) 
    elif gen == "uniform":
        #f = open("stat_r/" + str(n) +"/ntru_uniform/"+str(n) + "_" + str(beta), "w")
        filename =  "data-bkz1/"+"dimension" +str(n)+ "/" + "ntru_uniform/" +"r_i/"+str(times)+"times" +str(n)+ "n" +"_" + "q" + str(q) + "_"+ "beta"+ str(beta)  
    datamatrix = dat_to_matrix(filename, cols = n-1, divided_ch = '  ')
    #print(datamatrix[0] [139]  )          
    L1 = []
    L2 = []
    for m in range(len(datamatrix)):
        L1 += [float(datamatrix[m][i] )]
        L2 += [float(datamatrix[m][j] )]
        
    e1 = sum(L1) / len(L1)
    e2 = sum(L2) / len(L2)
    s = 0.0
    for i in range(len(L1)):
        s += (L1[i] - e1) * (L2[i] - e2)
    s /= len(L1)
    return s





def gen_stat(q,l_n, l_beta,gen,times):
    """Return the statistics of r_i's."""
    for n in l_n:
        for beta in l_beta:
            if gen == "gen_ntru_instance_circulant":
                f = open("stat_r/" + str(n) +"/ntru_circulant/"+str(n) + "_" + str(beta), "w")
                filename = "data-bkz1/"+"dimension" +str(n)+ "/" + "ntru_circulant/" + "r_i/"+str(times)+"times" +str(n)+ "n" +"_" + "q" + str(q) + "_"+ "beta"+ str(beta)  
            elif gen == "gen_ntru_instance_matrix":
                f = open("stat_r/" + str(n) +"/ntru_matrix/"+str(n) + "_" + str(beta), "w")
                filename =  "data-bkz1/"+"dimension" +str(n)+ "/" + "ntru_matrix/" + "r_i/"+str(times)+"times" +str(n)+ "n" +"_" + "q" + str(q) + "_"+ "beta"+ str(beta) 
            elif gen == "uniform":
                f = open("stat_r/" + str(n) +"/ntru_uniform/"+str(n) + "_" + str(beta), "w")
                filename =  "data-bkz1/"+"dimension" +str(n)+ "/" + "ntru_uniform/" +"r_i/"+str(times)+"times" +str(n)+ "n" +"_" + "q" + str(q) + "_"+ "beta"+ str(beta)  
            datamatrix = dat_to_matrix(filename, cols = n-1, divided_ch = '  ')
            for j in range(0, n-1):
                L = []
                for i in range(0, n):
                    L += [(datamatrix[i][j])]
                #print(L)
                #time.sleep(5)
                av, std = stat(L)
                f.write(str(av) + "  " + str(std) + "\n")
            f.close()



q = 31
l_beta = range(5, 31, 5)
#l_beta =[20]
l_n = [140]
n = 140
#gen_stat(q,l_n, l_beta, gen= "uniform",times=5000 )

#s= gen_cov(n= 140, beta=5, i=70, j=72, q=31 ,gen ="uniform" ,times=5000 )
#print(s)
#s= gen_cov(n= 140, beta=5, i=70, j=71, q=31 ,gen ="uniform" ,times=5000 )
#print(s)


# Plot Matrix
Beta = [20]

n = 140
Matrix = [[] for i in range(n )]


for beta in Beta:
    filename = "/home/zhaozishen/zhao-Z-shape/cov_mat_140_20"
    datamatrix = dat_to_matrix(filename, cols = n-1, divided_ch = ' ') 
    print( len(datamatrix) )
    print(len(datamatrix[0]  )  )
    plt.figure(figsize=(8, 5))
    plt.title("Covariance Matrix for BKZ_" + str(beta) + " Basis" +"uniform", fontsize=20)
    #X = np.array(datamatrix)
    cmap = matplotlib.cm.gray_r
    plt.pcolor(datamatrix, cmap=cmap)# edgecolor='white')
    plt.colorbar()
    plt.savefig("no_diag_cov_mat_uniform" + str(n) + "_" + str(beta) + ".png")




'''
f = open('cov_mat_' +"gen_ntru_instance_circulant"+ str(n) + '_' + str(beta), 'w')
for row in range(0, n-1):
        Row = []
        for col in range(0, n-1):
            if col != row:
                Row += [str(gen_cov(n, beta, row, col, q=31 ,gen ="gen_ntru_instance_circulant" ,times=5000)) + str(' ')]
                Matrix[row - 1] += [gen_cov(n, beta, row, col, q=31 ,gen ="gen_ntru_instance_circulant" ,times=5000)]
            else:
                Row += [str(0) + str(' ')]
                Matrix[row - 1] += [0]
        f.writelines(Row)
        f.write('\n')
'''

'''
for beta in Beta:
    Matrix = [[] for i in range(n )]
    f = open('cov_mat_' + str(n) + '_' + str(beta), 'w')
    for row in range(0, n-1):
        Row = []
        for col in range(0, n-1):
            if col != row:
                Row += [str(gen_cov(n, beta, row, col, q=31 ,gen ="gen_ntru_instance_circulant" ,times=5000)) + str(' ')]
                Matrix[row - 1] += [gen_cov(n, beta, row, col, q=31 ,gen ="uniform" ,times=5000)]
            else:
                Row += [str(0) + str(' ')]
                Matrix[row - 1] += [0]
        f.writelines(Row)
        f.write('\n')
    plt.figure(figsize=(8, 5))
    plt.title("Covariance Matrix for BKZ_" + str(beta) + " Basis", fontsize=20)
    X = np.array(Matrix)
    cmap = matplotlib.cm.gray_r
    plt.pcolor(X, cmap=cmap,edgecolor='white')
    plt.colorbar()
    plt.savefig("no_diag_cov_mat_" + str(n) + "_" + str(beta) + ".png")
'''

'''
no problem
for beta in Beta:
    filename = 'cov_mat_' + str(n) + '_' + str(beta)
    datamatrix = dat_to_matrix(filename, cols = n-1, divided_ch = ' ') 
    plt.figure(figsize=(8, 5))
    plt.title("Covariance Matrix for BKZ_" + str(beta) + " Basis", fontsize=20)
    X = np.array(datamatrix)
    cmap = matplotlib.cm.gray_r
    plt.pcolor(X, cmap=cmap, edgecolor='white')
    plt.colorbar()
    plt.savefig("no_diag_cov_mat_" + str(n) + "_" + str(beta) + ".png")
'''