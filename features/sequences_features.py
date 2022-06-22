def mutation_aa_label(three_letter):
 '''Label the residue'''
 aa = {'A':'1','R':'2','N':'3','D':'4','C':'5','E':'6','Q':'7','G':'8','H':'9','I':'10','L':'11','K':'12','M':'13','F':'14','P':'15','S':'16','T':'17','W':'18','Y':'19','V':'20'}
 return (aa.get(three_letter,0))

def sequences_features(chainA,mutation_resid,wild_aa, dir):
 '''This function labels the neighbourhood of the mutation site. It labels five residue behind the mutation position,
 mutation position itself, and five residue after the mutation point. This function is aquired from SAAFEC-SEQ,
 available at http://compbio.clemson.edu/SAAFEC-SEQ/download/SAAFEC-SEQ_code.zip
 Input: Chain identifier, mutant-residue, wild-residue, directory of the sequence file
 Output: Vector of 11 integers'''
 #4.2.3 Sequence Neighbours Feature
 mutation_resid=int(mutation_resid)-1
 features=[]
 target_sequence=""
 for line in open(dir+'/'+chainA):
  str=line.strip()
  if str[0]!=">":
   target_sequence=target_sequence+str
 #print(target_sequence)
 if len(target_sequence)<20:
  print("The fasta sequence length is at least 20!")
 for i in target_sequence:
  if i.upper() not in "ARNDCQEGHILKMFPSTWYV":
   print("Please check the input sequence!")
   sys.exit()
 if wild_aa.upper() != target_sequence[mutation_resid]:
  print("Wild type is not same as input sequence!")
  print (target_sequence[mutation_resid], mutation_resid)
  print (wild_aa)
  sys.exit()
 for i in range(int(mutation_resid)-5,int(mutation_resid)+6):
  if i < 0:
   features.append("0")
   continue
  try:
   target_sequence[i]
  except IndexError:
   features.append("0")
  else:
   features.append(mutation_aa_label(target_sequence[i]))
 #print('seq features#', len(features),features)
 return(features)