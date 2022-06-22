import numpy as np
#SPD3-numpy feature from BoostDDG
def parse_spd33(file_name, mutpos, deleted, dir):
    file_name += '.spd33'
    rnam1_std = "ACDEFGHIKLMNPQRSTVWY"
    ASA_std = (115, 135, 150, 190, 210, 75, 195, 175, 200, 170,
               185, 160, 145, 180, 225, 115, 140, 155, 255, 230)
    dict_rnam1_ASA = dict(zip(rnam1_std, ASA_std))
    with open(dir+'/'+file_name, 'r') as opf:
        lines=opf.read().split('\n')
        mut_line=lines[mutpos].rstrip() #pos is 1-baesd
        mut_line_sp=mut_line.split()
        #Struct position used as a feature
        #Struct=int(SS.get(mut_line_sp[2]))
        line = list(map(float, mut_line_sp[3:13]))
        temp_ASA_std = dict_rnam1_ASA[deleted]
        ASA = [line[0]/temp_ASA_std]
        HCEprob = np.array(line[7:10])
        ASA_SS = np.concatenate([ASA, HCEprob])
        ASA_SS = list(ASA_SS)
    return ASA_SS