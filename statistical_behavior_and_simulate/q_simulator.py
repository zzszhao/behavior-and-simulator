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




def gso (d , m , q , scale , prec =53):
    """
    Initial GSO of the q- ary lattice
    : param d: Lattice dimension . : param m: Number of ‘q‘ vectors
    : param q: Modulus ‘q‘
    """
    RR = RealField ( prec )
    r = [0] * d
    for i in range (0 , m , 1):
        r [ i ] = RR ( q ** 2)
    for j in range (m , d , 1):
        r [ j ] = scale ** 2
    return r
@cached_function
def ah_constant ( beta , scale , prec =53):
    RR = RealField ( prec )
    h = [
        RR (1.00000000000000) ,
        RR (1.00000000000000) ,
        RR (1.02492695474985) ,
        RR (1.05161980747554) ,
        RR (1.10445074764426) ,
        RR (1.12738094965280) ,
        RR (1.18113914825808) ,
        RR (1.21037399094975) ,
        RR (1.25380229130682) ,
        RR (1.31275912013118) ,
        RR (1.32737148909975) ,
        RR (1.36036038598451) ,
        RR (1.43966294868932) ,
        RR (1.46160335381056) ,
        RR (1.55544109937333) ,
        RR (1.51556419271870) ,
        RR (1.62261393856106) ,
        RR (1.62261542586149) ,
        RR (1.69386128187285) ,
        RR (1.74294139717984) ,
        RR (1.76511645842542) ,
        RR (1.85062593391314) ,
        RR (1.92910191727097) ,
        RR (1.95924301565493) ,
        RR (2.01827070103134) ,
        RR (2.09555655837831) ,
        RR (2.10779553489076) ,
        RR (2.16072906694966) ,
        RR (2.17866707693959) ,
        RR (2.24637113675736) ,
        RR (2.30374935080714) ,
        RR (2.32703241781875) ,
        RR (2.38719370287110) ,
        RR (2.44070021621789) ,
        RR (2.47587251254528) ,
        RR (2.57636556347475) ,
        RR (2.59580776592391) ,
        RR (2.70737106887894) ,
        RR (2.71293172456551) ,
        RR (2.80972868887803) ,
        RR (2.85337256682321) ,
        RR (2.96643112072014) ,
        RR (2.92901889083372) ,
        RR (2.99441283605227) ,
        RR (3.03884800785631) ,
        RR (3.09367293382578) ,
        RR (3.15593321975939) ,
        RR (3.21699662499604) ,
        RR (3.26306485591489) ,
        RR (3.32320705128515) ,
        RR (3.37871256503968) ,
        RR (3.49314144958915) ,
        RR (3.52315959119194) ,
        RR (3.60596719185516) ,
        RR (3.62302348800297) ,
        RR (3.67106448779994) ,
        RR (3.75288918941554) ,
        RR (3.80729231066355) ,
        RR (3.83017943995872) ,
        RR (3.87680907526471) ,
        RR (3.99304695483873) ,
        RR (4.04539449313441) ,
        RR (4.12835783024051) ,
        RR (4.12989317768338) ,
        RR (4.25784341671535) ,
        RR (4.23117443161792) ,
        RR (4.34900834929189) ,
        RR (4.37617075735239) ,
        RR (4.44178322064643) ,
        RR (4.53471153665029) ,
        RR (4.54352770842019) ,
        RR (4.65495019943735) ,
        RR (4.63508961564190) ,
        RR (4.73512141747094) ,
        RR (4.76413486921241) ,
        RR (4.87996877805749) ,
        RR (4.94652066938356) ,
        RR (4.98495063065338) ,
        RR (5.07692160777643) ,
        RR (5.04745948303442) ,
        RR (5.12518339690435) ,
        ]


    ah = [ RR (0.0)] * ( beta + 1)
    c = [ RR (0.0)] * ( beta + 1)
    for j in range (0 , min ( beta , 80) + 1 , 1):
        ah [ j ] = h [ j ]

    
    if beta > 49:
        lga = RR (
            log_gamma ( beta / 2 + 1) * (2 / beta ) - log ( pi )
        )
        ah [ beta ] = RR ( exp ( lga ))
        # ah[ beta ]= GH( beta )ˆ2. x is chosen to smooth ah[j] for j in [x,beta -1] ,
        # where ’ ’7 ’ ’ is the ’’best ’’ using the above data h[i] ’s.
        z = min (
            7 + floor ( log ( scale , 2)) , beta - 4
        ) # z is related to scale , but can be upper bounded .
        x = min ( beta - z , 80)
        ah [ x - 1] = (
            ( ah [ beta ] - h [ x - 2]) / ( beta - x + 2) / 2
            + h [ x - 2] / 2
            + h [ x - 1] / 2
        )
        for j in range (x , beta , 1):
            # (j,u) is on the stright line through two points (x -2 ,h[x -2]) and (beta ,ah[ beta ]).
            u = ( ah [ beta ] - h [ x - 2]) * ( j - x + 2) / ( beta - x + 2) + h [ x - 2]
            # (j,v) is on the stright line through two points (x -1 ,h[x -1]) and (beta ,ah[ beta ]).
            v = ( ah [ beta ] - h [ x - 1]) * ( j - x + 1) / ( beta - x + 1) + h [ x - 1]
            # (j,w) is on the straight line through two points (x,h[x]) and (beta ,ah[ beta ]).
            w = ( ah [ beta ] - h [ x ]) * ( j - x ) / ( beta - x ) + h [ x ]
            ah [ j ] = (
                u + v + w
            ) / 3 # We take the average of u, v and w as ah[j].
    c [0] = 0
    c [1] = 0
    for j in range (2 , beta + 1 , 1):
        c [ j ] = 0.5 * RR ( log ( ah [ j ] , 2.0))  # log (sqrt(r))
    return c





def one_tour (d, m, q, r , block_size , c ):
    """
    Updates the log GSO for one ( simulated ) BKZ tour on q- ary lattice .
    """
    d = len ( r )
    r1 = copy(r)
    tail_length = int ( ( m - int( 1/2 + ln(q) / (2/(block_size -1)* (c[block_size]/log(2,e))) )  ) /4 )
    X = random.expovariate(.5)
    
    #print( tail_length)
    '''
    if tail_length >0 :
        block = 0
        logX = 0
    else :
        block = 45
        logX = log(X,2)
    '''
    for i in range ( d - 1  ):
        p = min( block_size , d - i ) # beta_i
        j = min( i + block_size , d ) # n_i
        g = sum ([ r [ z ] for z in range (i , j , 1)]) # alogrithm5 line13
        if r [ i ] > c [ p ] + ( g   ) / p   :  # bsw_c[p-1]  replace  c[p]
            r [ i ] =c [ p ] + ( g    ) / p
            for s in range ( i + 1 , j , 1):
                t = j - s
                y = sum ([ r [ z ] for z in range (i , s , 1)])
                r [ s ] = c [ t ] + ( g - y   ) / t
    return r


def qary_simulate_r (d, m, q, r , block_size , tours =1 , scale =1 , prec =53):
    RR = RealField ( prec )
    c = ah_constant ( block_size , scale )
    r = [ RR ( log ( r_ , 2.0)) / 2.0 for r_ in r ]
    for i in range ( tours ):
        r = one_tour (d, m, q, r , block_size , c )
        #file =   open( "20220826-22:46RHF-smiluator"+str(tours)+str("-unifrom")  +str(d)+"-"+str()+"-"+str(q)+"-" +"beta"+str(beta) ,'a' )
        logRHF = ( r[0] - (log(q,2)/2.0 )   )/d
        RHF = 2**logRHF
        #print(float(RHF))
        #file.writelines(RHF)
        #file.write('\n')
        #file.close() 
    r = [2.0 ** (2 * r_ ) for r_ in r ]   # |b_j^{*}| ^{2}
    return r


def qary_simulate (d , m , q , block_size , tours =1 , scale =1):
    r = gso (d , m , q , scale )
    r = qary_simulate_r (d, m, q, r , block_size = block_size , tours = tours , scale = scale )
    
    return r



def normalize_GSO_unitary(l):
    log_det = sum(l)
    n = len(l)
    nor_log_det = [0.0] * n
    for i in range(n):
        nor_log_det[i] = l[i] - log_det/n
    return nor_log_det



def zsimulate_prob(d, m, q, r, param,  c, prng_seed=88 ):

    if param.block_size <= 2:
        raise ValueError("The simulator requires block size >= 3.")
    # fix PRNG seed
    random.seed(prng_seed if prng_seed else FPLLL.randint(0, 2**32-1))

    #r = _extract_log_norms(r)
    r1 = copy(r)
    r2 = copy(r)

    if param.max_loops:
        N = param.max_loops
    else:
        N = d
    beta = param.block_size
    q_length = int ( ( m - int( 1/2 + ln(q) / (2/(beta-1)* (c[beta]/log(2,e))) )  ) /4 )
    #print(q_length)
    t0 = [True for _ in range(d)]
    for i in range(N):
        t1 = [False for _ in range(d)]
        #first 
         #first 
        for k in range(d - min(50, param.block_size)):  #for k in range(d - min(45, param.block_size)):
            beta = min(param.block_size, d - k)
            f = k + beta
            phi = False
            for kp in range(k, f):
                phi |= t0[kp]
            logV = sum(r1[:f]) - sum(r2[:k])
            if phi:
                X = random.expovariate(.5)
                lma = (log(X, 2) + logV) / beta + c[beta]
                if lma < r1[k]:
                    r2[k] = lma
                    r2[k+1] = r1[k] + log(sqrt(1-1./beta), 2)
                    dec = (r1[k]-lma) + (r1[k+1] - r2[k+1])
                    for j in range(k+2, f):
                        r2[j] = r1[j] + dec/(beta-2.)
                        t1[j] = True
                    phi = False

            for j in range(k, f):
                r1[j] = r2[j]
        # early termination
        if True not in t1:
            break
       

        beta = min( 50 ,  param.block_size)
        logV = sum(r1) - sum(r2[0:d-beta])


        for k in range(d-beta, d):
            r2[k] = logV / (d-k) + c[d-k]
            logV = logV - r2[k]
            t1[k] = True
        logV = sum(r) -sum(r2)
        #print(logV)
        #print("the final phase")

        if (r1 == r2):
            break
        r1 = copy(r2)
        t0 = copy(t1)
        #file =   open( "20220826-22:46RHF-ourzsmiluator"+str(tours)+str("-unifrom")  +str(N)+"-"+str()+"-"+str(q)+"-" +"beta"+str(beta) ,'a' )
        logRHF = ( r1[0] - (log(q,2)/2.0 )   )/d
        RHF = 2**logRHF
        print( float(RHF))
        #file.writelines(RHF)
        #file.write('\n')
        #file.close()      
    #else :
    return r1 




def pro_qary_simulate_r (d,m,q, r , block_size , tours =1 , scale =1 , prec =53):
    RR = RealField ( prec )
    c = ah_constant ( block_size , scale )
    r = [ RR ( log ( r_ , 2.0)) / 2.0 for r_ in r ]
    r1 = zsimulate_prob(d, m, q, r, BKZ.Param(block_size=block_size, max_loops=tours, flags=BKZ.VERBOSE) ,c   )
    r1 = [2.0 ** (2 * r_ ) for r_ in r1 ]   # |b_j^{*}| ^{2}

    return r1


def pro_qary_simulate (d , m , q , block_size , tours =1 , scale =1):
    r = gso (d , m , q , scale )
    r1 = pro_qary_simulate_r (d,m, q, r , block_size = block_size , tours = tours , scale = scale )
    return r1



