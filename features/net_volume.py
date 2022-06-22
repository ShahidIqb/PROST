def net_volume(wild,mutation):
 '''Subtract the wild-residue volume from mutant-residue volume. This function is aquired from SAAFEC-SEQ,
 available at http://compbio.clemson.edu/SAAFEC-SEQ/download/SAAFEC-SEQ_code.zip'''
 vol = {'A':'-0.594','R':'0.316','N':'-0.32','D':'-0.353','C':'-0.381','E':'-0.06','Q':'-0.002','G':'-0.9','H':'0.099',
        'I':'0.244','L':'0.244','K':'0.265','M':'0.203','F':'0.493','P':'-0.335','S':'-0.59','T':'-0.299','W':'0.9',
        'Y':'0.533','V':'-0.042'}
 return ('{:.1f}'.format(float(vol.get(mutation,'0'))-float(vol.get(wild,'0'))))