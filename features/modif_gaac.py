def modif_gaac(wild,mutation):
 '''This function gives a Grouped amino-acid composition label according to the wild-type and mutant residue.
 This function is acquired from iFeature, available at https://github.com/Superzchen/iFeature
 Input: wild-residue, mutant-residue
 Output: Integer label'''
 res_grp1 = ('G','A','V','L','M','I')
 res_grp2 = ('F','Y','W')
 res_grp3 = ('K','R','H')
 res_grp4 = ('D','E')
 res_grp5 = ('S','T','C','P','N','Q')
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