import numpy as np

def get_blosum_line(wt):
    aa_index = {'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4, 'Q': 5, 'E': 6, 'G': 7, 'H': 8, 'I': 9, 'L': 10, 'K': 11,
                'M': 12, 'F': 13, 'P': 14, 'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19}
    f=open('blosum62.iij','r')
    lines=f.read().split('\n')
    line_idx=aa_index[wt]+2
    line=lines[line_idx]
    line_sp=line.split()
    blosum_line=list(map(float,line_sp[1:21]))
    blosum_array=np.array(blosum_line)
    return blosum_array

#dF score from BoostddG method, uses pssm score to calculate the difference betwween wild and mutant-type
def get_dF(file_name, mutpos, deleted, induced, dir):
    aa_index = {'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C': 4, 'Q': 5, 'E': 6, 'G': 7, 'H': 8, 'I': 9, 'L': 10, 'K': 11,
                'M': 12, 'F': 13, 'P': 14, 'S': 15, 'T': 16, 'W': 17, 'Y': 18, 'V': 19}
    file_name += '.pssm'
    with open(dir+'/'+file_name, 'r') as opf:
        lines=opf.read().split('\n')
        #print(file_name, mutpos, deleted, induced)
        line_idx=mutpos+2 #pos is 1-based
        line=lines[line_idx].strip()
        line_sp=line.split()
        F_line=list(map(float, line_sp[22:42]))
        F_array=np.array(F_line)
        if np.all(F_array==0):
            F_array=get_blosum_line(deleted)
            F_array=F_array/10
        else:
            F_array=F_array/100
        Fwt=F_array[aa_index[deleted]]
        Fmt=F_array[aa_index[induced]]
        dF=Fmt-Fwt
        return float(dF)