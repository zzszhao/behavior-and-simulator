"""
The python code is adapted from yang yu and Léo Ducas's implementation.
"""

## generate three lattice types ::: ntru_matrix, ntru_circulant, q-ary-uniform
## record the result:  L_log  L_rate_bi_GH   L_r_log

# L_log :::  the log of  |b_j^{*}| - average_det_log
# L_rate_bi_GH  :: the log of  |b_j^{*}| / vol(B[j, min(j+beta, n)])
# L_r_log   :::   log value  |b_j^{*}| / |b_{j+1}^{*}|


from ntru_keygen import gen_ntru_instance_matrix, gen_ntru_instance_circulant
from fpylll import IntegerMatrix, BKZ, FPLLL, LLL, GSO
from numpy import array, zeros, identity, block

from math import log, sqrt, lgamma, pi, exp
import numpy as np

from fpylll.algorithms.bkz2 import BKZReduction as BKZ2
from fpylll.algorithms.bkz import BKZReduction as BKZ1

from fpylll.algorithms.simple_bkz import BKZReduction as simpleBKZ

import time


sphere = [        (lgamma(beta / 2.0 + 1) * (1.0 / beta) - log(sqrt(pi))) / log(2.0)
        for beta in range(1, 200 + 1)     ]   #GH(beta)  

#  list(map(lambda x: log(x, 2) / 2.0, d  ))

def generate_data(n, q, betasize, times, bkz, gen ,seed , tours ):
    N = 2*n
    average_det_log = n*log(q , 2 )/N
    sigmasq = 0.667
    Block = range(5, betasize+1, 5) #output results when beta=
    Block_progressive = range(2, betasize+1, 1 )
    B = IntegerMatrix(N, N )
    if gen == "gen_ntru_instance_circulant":
        file_list =  [ open("data-bkz1/" +"dimension" +str(N)+ "/"+ "ntru_circulant/" + "square_norm/"+str(times)+"times"+str(N) +"n" +"_" + "q" + str(q) + "_"+ "beta"+ str(beta) ,'a' )  for beta in Block ]
        file_r_list =  [ open("data-bkz1/"+"dimension" +str(N)+ "/" + "ntru_circulant/" + "r_i/"+str(times)+"times" +str(N)+ "n" +"_" + "q" + str(q) + "_"+ "beta"+ str(beta) ,'a' )  for beta in Block ]
        file_rate_list =  [ open("data-bkz1/"+"dimension" +str(N)+ "/" + "ntru_circulant/" + "rate_bi_GH/"+str(times)+"times" +str(N)+ "n" +"_" + "q" + str(q) + "_"+ "beta"+ str(beta) ,'a' )  for beta in Block ]
        print("oepn circulant ntru")
    elif gen == "gen_ntru_instance_matrix":
        file_list =  [ open("data-bkz1/"+"dimension" +str(N)+ "/" + "ntru_matrix/" + "square_norm/"+str(times)+"times" +str(N)+ "n" +"_" + "q" + str(q) + "_"+ "beta"+ str(beta) ,'a' )  for beta in Block ]
        file_r_list =  [ open("data-bkz1/"+"dimension" +str(N)+ "/" + "ntru_matrix/" + "r_i/"+str(times)+"times" +str(N)+ "n" +"_" + "q" + str(q) + "_"+ "beta"+ str(beta) ,'a' )  for beta in Block ]
        file_rate_list =  [ open("data-bkz1/"+"dimension" +str(N)+ "/" + "ntru_matrix/" + "rate_bi_GH/"+str(times)+"times" +str(N)+ "n" +"_" + "q" + str(q) + "_"+ "beta"+ str(beta) ,'a' )  for beta in Block ]
        print("open matrix ntru")
    elif gen == "uniform":
        file_list =  [ open("data-bkz1/"+"dimension" +str(N)+ "/" + "ntru_uniform/" +"square_norm/"+str(times)+"times" +str(N)+ "n" +"_" + "q" + str(q) + "_"+ "beta"+ str(beta) ,'a' )  for beta in Block ]
        file_r_list =  [ open("data-bkz1/"+"dimension" +str(N)+ "/" + "ntru_uniform/" +"r_i/"+str(times)+"times" +str(N)+ "n" +"_" + "q" + str(q) + "_"+ "beta"+ str(beta) ,'a' )  for beta in Block ]
        file_rate_list =  [ open("data-bkz1/"+"dimension" +str(N)+ "/" + "ntru_uniform/" +"rate_bi_GH/"+str(times)+"times" +str(N)+ "n" +"_" + "q" + str(q) + "_"+ "beta"+ str(beta) ,'a' )  for beta in Block ]
        print("open unifrom")
    for s in range(times):
        if gen == "gen_ntru_instance_matrix":
            A, F, G = gen_ntru_instance_matrix(n, q, sigmasq, seed+s)
        elif gen == "gen_ntru_instance_circulant":
            A, F, G = gen_ntru_instance_circulant(n, q, sigmasq, seed+s)
        elif gen == "uniform":
            D = IntegerMatrix(n, n )
            FPLLL.set_random_seed( seed+s )
            D.randomize("uniform", bits=5)
            DD = np.zeros((n, n))
            for i in range(n):
                for j in range(n):
                    DD[i][j] = D[i,j]
            A = block([[q * identity(n, dtype="long") , zeros((n, n), dtype="long") ],      [   DD  , identity(n, dtype="long") ] ])
        for i in range(N):
            for j in range(N):
                B[i,j] = int(A[i][j])
        # print(B)
        bkz_B = simpleBKZ(B)     
        for i in range(len(Block_progressive)):
            beta = Block_progressive[i]
            #param = BKZ.Param(block_size = beta, max_loops=tours, strategies = BKZ.DEFAULT_STRATEGY)
            if bkz == "bkz1":
                bkz_B.__call__(beta)
                param = BKZ.Param(block_size = beta, max_loops=tours, strategies = BKZ.DEFAULT_STRATEGY)
                #BKZ1(B)( param )
                # print("bkz1::::::::")
            else :
                param = BKZ.Param(block_size = beta, max_loops=tours, strategies = BKZ.DEFAULT_STRATEGY)
                BKZ2(B)( param )
            gso_B = GSO.Mat(B)
            gso_B.update_gso()
            #L = [str(gso_B.get_r(j, j) ) + "  " for j in range(N)]
            L_log = [str( log(gso_B.get_r(j, j), 2)/2.0 - average_det_log   ) + "  " for j in range(N)]
            # Ignore the infulence of lattice volume  by subtract average_det_log
            #
            if beta in Block:  
            #file_list[i].writelines(L)
                index = int(beta/5)-1
                file_list[index].writelines(L_log)
                file_list[index].write("\n")
                #L_log = [str( log(gso_B.get_r(j, j), 2)/2.0 - average_det_log   ) + "  " for j in range(N)]
                log_GH = [float(0)]*N
                for k in range(N):
                    log_vol_beta = 0.0
                    if beta >= N-k :
                        for kk in range(N-k):
                            log_vol_beta = log_vol_beta + float(L_log[N-kk-1])
                        log_GH[k] = (1.0 / (N-k)) * log_vol_beta + sphere[N-k-1]
                    else:
                        for kk in range(beta):
                            log_vol_beta = log_vol_beta + float(L_log[k+kk])
                        log_GH[k] = (1.0 / beta) * log_vol_beta + sphere[beta-1]

                L_rate_bi_GH = [ str(  float(L_log[j]) - log_GH[j]  ) + "  " for j in range(N) ]
                file_rate_list[index].writelines(L_rate_bi_GH)
                file_rate_list[index].write("\n")
                L_r_log = [ str(  float(L_log[j])  - float(L_log[j+1])  )+ "  " for j in range(N-1)]
                # (log_b_i* - average_det _log ) - (log_b_{i+1}* - average_det _log )
                file_r_list[index].writelines(L_r_log)
                file_r_list[index].write("\n")
          
#t = time.perf_counter()
#run(n, q , 0, 10, 1 ,bkz="bkz1", gen="uniform"  )
#print(" Time for (BKZ1) blocksize "  + " uniform" + str(n) + ":", time.perf_counter()-t, "s" )



#t = time.perf_counter()
#run(n, q , 0, 10, 1 ,bkz="bkz1", gen="gen_ntru_instance_matrix"  )
# print(" Time for (BKZ1) blocksize "  + "gen_ntru_instance_matrix" + str(n) + ":", time.perf_counter()-t, "s" )
