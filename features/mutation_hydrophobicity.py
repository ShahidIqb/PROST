def mutation_hydrophobicity(wild,mutation):
 '''This function returns a hydrophobicity label according to the wild-type and mutant residue. This function is aquired
 from SAAFEC-SEQ, available at http://compbio.clemson.edu/SAAFEC-SEQ/download/SAAFEC-SEQ_code.zip
 Input: wild-residue, mutant-residue
 Output: Integer label'''
 res_grp1 = ('I','V','L','F','C','M','A','W')
 res_grp2 = ('G','T','S','Y','P','H')
 res_grp3 = ('N','D','Q','E','K','R')
 if wild in res_grp1:
  if mutation in res_grp2:
   return 0
  elif mutation in res_grp3:
   return 1
  elif mutation in res_grp1:
   return 2
 elif wild in res_grp2:
  if mutation in res_grp2:
   return 3
  elif mutation in res_grp3:
   return 4
  elif mutation in res_grp1:
   return 5
 elif wild in res_grp3:
  if mutation in res_grp2:
   return 6
  elif mutation in res_grp3:
   return 7
  elif mutation in res_grp1:
   return 8