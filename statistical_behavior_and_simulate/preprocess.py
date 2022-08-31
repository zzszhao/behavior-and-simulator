# -* - coding : utf -8 -* -
# 

from sage . all import e,log ,ln, log_gamma , pi , exp , floor , RealField

from sage . all import cached_function
from copy import copy
import numpy as np
import matplotlib.pyplot as plt

from  read_data  import dat_to_one_list
from fpylll import IntegerMatrix, BKZ, FPLLL, LLL, GSO
from numpy import array, zeros, identity, block

from math import sqrt, lgamma, pi, exp
import numpy as np

from fpylll.algorithms.bkz2 import BKZReduction as BKZ2
from fpylll.algorithms.bkz import BKZReduction as BKZ1
from ntru_keygen import gen_ntru_instance_matrix, gen_ntru_instance_circulant
import time, random 







#average_det_log = 80*log(17 , 2 )/180

#print("***************")
#print(sum(pro_r_log_det1))

#print( sum(r_log_det1) - sum(pro_r_log_det1)  )
#time.sleep(5)


def ceshi(d,m,q,gen,beta,tours,seed=15, sigmasq = float(2/3) ):
    B = IntegerMatrix(d, d )
    average_det_log = m*log(q , 2 )/d
    if gen=="uniform":
        D = IntegerMatrix(d, d )
        FPLLL.set_random_seed( seed )
        D.randomize("uniform", bits=5)
        DD = np.zeros((d-m, m))
        for i in range(d-m):
            for j in range(m):
                DD[i][j] = D[i,j]
        A = block([[q * identity(m, dtype="long") , zeros((m, d-m), dtype="long") ],      [   DD  , identity(d-m, dtype="long") ] ])
    if gen == "circulant":
        A, F, G = gen_ntru_instance_circulant(m, q, sigmasq, seed)
    for i in range(d):
        for j in range(d):
            B[i,j] = int(A[i][j])
    for betasize in range(3,beta):
        param = BKZ.Param(block_size = betasize, max_loops=tours, strategies = BKZ.DEFAULT_STRATEGY)
        BKZ2(B)( param )
    gso_B = GSO.Mat(B)
    gso_B.update_gso()
    B_log_det1 = [float( log(gso_B.get_r(j, j), 2)/2.0 - average_det_log   )  for j in range(d)]
    return B_log_det1



