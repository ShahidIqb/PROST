import numpy as np
import math

def cal_norm_aap(row,col):
    aap_list = [[0.5,-0.5,89.3,115,0.305,1.29,1.08],
                [0,3,190.3,225,0.227,1,1.05],
                [0,0.2,122.4,160,0.322,0.81,0.85],
                [0,3,114.4,150,0.335,1.1,0.85],
                [0,-1,102.5,135,0.339,0.79,0.95],
                [0,0.2,146.9,180,0.306,1.07,0.95],
                [0,3,138.8,190,0.282,1.49,1.15],
                [0,0,63.8,75,0.352,0.63,0.55],
                [0.5,-0.5,157.5,195,0.215,1.33,1],
                [1.8,-1.8,163,175,0.278,1.05,1.05],
                [1.8,-1.8,163.1,170,0.262,1.31,1.25],
                [0,3,165.1,200,0.391,1.33,1.15],
                [1.3,-1.3,165.8,185,0.28,1.54,1.15],
                [2.5,-2.5,190.8,210,0.195,1.13,1.1],
                [0,0,121.6,145,0.346,0.63,0.71],
                [0,0.3,94.2,115,0.326,0.78,0.75],
                [0.4,-0.4,119.6,140,0.251,0.77,0.75],
                [3.4,-3.4,226.4,255,0.291,1.18,1.1],
                [2.3,-2.3,194.6,230,0.293,0.71,1.1],
                [1.5,-1.5,138.2,155,0.291,0.81,0.95]]
    aap_array = np.array(aap_list,dtype=float)
    avg_col_sum = np.sum(aap_array,axis=0)/20
    Iji_hat = aap_array[row][col]
    avg_Ij_hat = avg_col_sum[col]
    sumI = 0
    for k in range(20):
        Ijk_hat = aap_array[k][col]
        dI = (Ijk_hat - avg_Ij_hat)*(Ijk_hat - avg_Ij_hat)
        sumI += dI
    sumI = sumI/20
    root_sumI = math.sqrt(sumI)
    Iji = (Iji_hat-avg_Ij_hat)/root_sumI
    return Iji
#aa_score from BoostddG
def cal_aa_score(deleted,induced):
    aa_index = {'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4, 'Q': 5, 'E': 6, 'G': 7, 'H': 8, 'I': 9, 'L': 10, 'K': 11,
                'M': 12, 'F': 13, 'P': 14, 'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19}
    aa1_row = aa_index[deleted]
    aa2_row = aa_index[induced]
    sumI = 0
    for j in range(7):
        Ijaa1 = cal_norm_aap(aa1_row,j)
        Ijaa2 = cal_norm_aap(aa2_row,j)
        dIj = (Ijaa1-Ijaa2)*(Ijaa1-Ijaa2)
        sumI += dIj
    return sumI/7