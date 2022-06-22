def mutation_size(wild,mutation):
 '''This function returns a size label according to the wild-type and mutant residue. This function is aquired
 from SAAFEC-SEQ, available at http://compbio.clemson.edu/SAAFEC-SEQ/download/SAAFEC-SEQ_code.zip
 Input: wild-residue, mutant-residue
 Output: Integer label'''
 res_grp1 = ('G','A','S')
 res_grp2 = ('C','D','P','N','T')
 res_grp3 = ('E','V','Q','H')
 res_grp4 = ('M','I','L','K','R')
 res_grp5 = ('F','Y','W')
 if wild in res_grp1:
  if mutation in res_grp2:
   return 0
  elif mutation in res_grp3:
   return 1
  elif mutation in res_grp4:
   return 2
  elif mutation in res_grp5:
   return 3
  elif mutation in res_grp1:
   return 4
 elif wild in res_grp2:
  if mutation in res_grp2:
   return 5
  elif mutation in res_grp3:
   return 6
  elif mutation in res_grp4:
   return 7
  elif mutation in res_grp5:
   return 8
  elif mutation in res_grp1:
   return 9
 elif wild in res_grp3:
  if mutation in res_grp2:
   return 10
  elif mutation in res_grp3:
   return 11
  elif mutation in res_grp4:
   return 12
  elif mutation in res_grp5:
   return 13
  elif mutation in res_grp1:
   return 14
 elif wild in res_grp4:
  if mutation in res_grp2:
   return 15
  elif mutation in res_grp3:
   return 16
  elif mutation in res_grp4:
   return 17
  elif mutation in res_grp5:
   return 18
  elif mutation in res_grp1:
   return 19
 elif wild in res_grp5:
  if mutation in res_grp2:
   return 20
  elif mutation in res_grp3:
   return 21
  elif mutation in res_grp4:
   return 22
  elif mutation in res_grp5:
   return 23
  elif mutation in res_grp1:
   return 24