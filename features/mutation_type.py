def mutation_type(wild,mutation):
 '''This function labels wild-residue in accordance with the mutant-residue. This function is acquired from SAAFEC-SEQ,
 available at http://compbio.clemson.edu/SAAFEC-SEQ/download/SAAFEC-SEQ_code.zip
 Input: wild-residue, mutant-residue
 Output: Integer label'''
 wild_lists = ['A','F','C','D','N','E','Q','G','H','L','I','K','M','P','R','S','T','V','W','Y']
 mutation_lists = ['A','F','C','D','N','E','Q','G','H','L','I','K','M','P','R','S','T','V','W','Y']
 label_1=0
 label_2=0
 for i in wild_lists:
  for j in mutation_lists:
   if i != j:
    label_1 += 1
    if wild == i:
     if mutation == j and i != j:
      label_2 = label_1
      break
 return label_2