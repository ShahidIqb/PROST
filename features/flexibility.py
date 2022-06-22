def flexibility(wild,mutation):
 '''Subtract the flexibility scores of wild-residue from mutant-residue. This function is aquired from SAAFEC-SEQ,
 available at http://compbio.clemson.edu/SAAFEC-SEQ/download/SAAFEC-SEQ_code.zip'''
 flex = {'A':'-0.9','R':'0.446','N':'-0.311','D':'-0.614','C':'-0.866','E':'-0.008','Q':'0.9','G':'-0.9','H':'-0.311',
         'I':'-0.765','L':'-0.765','K':'0.446','M':'-0.463','F':'-0.614','P':'-0.883','S':'-0.866','T':'-0.866',
         'W':'-0.311','Y':'-0.614','V':'-0.866'}
 return (float(flex.get(mutation,'0'))-float(flex.get(wild,'0')))