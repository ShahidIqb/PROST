def mutation_chemical(wild,mutation):
 '''This function returns a label (chemical) according to the wild-type and mutant residue. This function is aquired
 from SAAFEC-SEQ, available at http://compbio.clemson.edu/SAAFEC-SEQ/download/SAAFEC-SEQ_code.zip
 Input: wild-residue, mutant-residue
 Output: Integer label'''
 res_grp1 = ('A','G','I','L','P','V')
 res_grp2 = ('R','H','K')
 res_grp3 = ('N','Q')
 res_grp4 = ('D','E')
 res_grp5 = ('C','M')
 res_grp6 = ('S','T')
 res_grp7 = ('F','W','Y')

 if wild in res_grp1:
  if mutation in res_grp2:
   return 0
  elif mutation in res_grp3:
   return 1
  elif mutation in res_grp4:
   return 2
  elif mutation in res_grp5:
   return 3
  elif mutation in res_grp6:
   return 4
  elif mutation in res_grp7:
   return 5
  elif mutation in res_grp1:
   return 6
 elif wild in res_grp2:
  if mutation in res_grp2:
   return 7
  elif mutation in res_grp3:
   return 8
  elif mutation in res_grp4:
   return 9
  elif mutation in res_grp5:
   return 10
  if mutation in res_grp6:
   return 11
  elif mutation in res_grp7:
   return 12
  elif mutation in res_grp1:
   return 13
 elif wild in res_grp3:
  if mutation in res_grp2:
   return 14
  elif mutation in res_grp3:
   return 15
  elif mutation in res_grp4:
   return 16
  elif mutation in res_grp5:
   return 17
  elif mutation in res_grp6:
   return 18
  elif mutation in res_grp7:
   return 19
  elif mutation in res_grp1:
   return 20
 elif wild in res_grp4:
  if mutation in res_grp2:
   return 21
  elif mutation in res_grp3:
   return 22
  elif mutation in res_grp4:
   return 23
  elif mutation in res_grp5:
   return 24
  elif mutation in res_grp6:
   return 25
  elif mutation in res_grp7:
   return 26
  elif mutation in res_grp1:
   return 27
 elif wild in res_grp5:
  if mutation in res_grp2:
   return 28
  elif mutation in res_grp3:
   return 29
  elif mutation in res_grp4:
   return 30
  elif mutation in res_grp5:
   return 31
  elif mutation in res_grp6:
   return 32
  elif mutation in res_grp7:
   return 33
  elif mutation in res_grp1:
   return 34
 elif wild in res_grp6:
  if mutation in res_grp2:
   return 35
  elif mutation in res_grp3:
   return 36
  elif mutation in res_grp4:
   return 37
  elif mutation in res_grp5:
   return 38
  elif mutation in res_grp6:
   return 39
  elif mutation in res_grp7:
   return 40
  elif mutation in res_grp1:
   return 41
 elif wild in res_grp7:
  if mutation in res_grp2:
   return 42
  elif mutation in res_grp3:
   return 43
  elif mutation in res_grp4:
   return 44
  elif mutation in res_grp5:
   return 45
  elif mutation in res_grp6:
   return 46
  elif mutation in res_grp7:
   return 47
  elif mutation in res_grp1:
   return 48