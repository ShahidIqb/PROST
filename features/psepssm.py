def psepssm(file_name, dir):
 import math
 lamda=7
 pssm=[]
 features=[]
 file_name += '.pssm'
 for line1 in open(dir+'/'+file_name, 'r'):
  info1=line1.strip().split()
  if len(info1)>43:
   if info1[1].isupper():
    length=int(info1[0])
    for i in range(2,22):
     pssm.append(float(1)/(1+math.e**(-int(info1[i]))))
 for i in range(20):
  sum_pssm=0
  for j in range(int(length)):
   sum_pssm=sum_pssm+float(pssm[i+20*j])
  features.append("%.2f"%(float(sum_pssm)/length))
 for i in range(1,lamda+1):
   for j in range(20):
    s_pssm=0
    for k in range(int(length)-i):
     s_pssm=(float(pssm[20*k+j])-float(pssm[20*(k+i)+j]))**2+s_pssm
    features.append("%.2f"%(float(s_pssm)/(int(length)-i)))
 #print('psePSSM#', len(features),features)
 return(features)