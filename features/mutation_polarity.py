def mutation_polarity(wild,mutation):
 '''This function returns a polarity label according to the wild-type and mutant residue. This function is aquired
 from SAAFEC-SEQ, available at http://compbio.clemson.edu/SAAFEC-SEQ/download/SAAFEC-SEQ_code.zip
 Input: wild-residue, mutant-residue
 Output: Integer label'''
 res_grp1 = ('R','H','K')
 res_grp2 = ('A','C','G','I','L','M','F','P','W','V')
 res_grp3 = ('N','Q','S','T','Y')
 res_grp4 = ('D','E')
 if wild in res_grp1:
  if mutation in res_grp2:
   return 0
  elif mutation in res_grp1:
   return 1
  elif mutation in res_grp3:
   return 2
  elif mutation in res_grp4:
   return 3
 elif wild in res_grp2:
  if mutation in res_grp2:
   return 4
  elif mutation in res_grp1:
   return 5
  elif mutation in res_grp3:
   return 6
  elif mutation in res_grp4:
   return 7
 elif wild in res_grp3:
  if mutation in res_grp2:
   return 8
  elif mutation in res_grp1:
   return 9
  elif mutation in res_grp3:
   return 10
  elif mutation in res_grp4:
   return 11
 elif wild in res_grp4:
  if mutation in res_grp2:
   return 12
  elif mutation in res_grp1:
   return 13
  elif mutation in res_grp3:
   return 14
  elif mutation in res_grp4:
   return 15