def net_hydrophobicity(wild,mutation):
 '''Subtract the hydrophobicity scores of wild-residue from mutant-residue. This function is aquired from SAAFEC-SEQ,
 available at http://compbio.clemson.edu/SAAFEC-SEQ/download/SAAFEC-SEQ_code.zip'''
 hyd = {'A':'-0.378','R':'0.502','N':'0.445','D':'0.321','C':'-0.262','E':'0.011','Q':'0.336','G':'0.03','H':'0.751',
        'I':'-0.748','L':'-0.808','K':'0.9','M':'-0.558','F':'-0.9','P':'-0.739','S':'0.056','T':'0.044','W':'-0.468',
        'Y':'-0.637','V':'-0.563'}
 return ('{:.1f}'.format(float(hyd.get(mutation,'0'))-float(hyd.get(wild,'0'))))