def mpssmscores(file_name,position, dir):
 #4.2.2 Neighbour Mutation Conservation Scores
 windows=7
 pssm=[]
 line1=[]
 file_name += '.pssm'
 for line2 in open(dir+'/'+file_name, 'r'):
  line1.append(line2.strip())
 for j in range(int(position)-(windows//2),int(position)+windows//2+1):
  if j<0:
   for i in range(2,22):
    pssm.append(0)
   continue
  index=0
  for info1 in line1:
   info2=info1.split()
   if len(info2)>43:
    if int(info2[0])==j:
     for i in range(2,22):
      pssm.append(info2[i])
     index=1
    if info1==line1[-7] and index==0:
     for i in range(2,22):
      pssm.append(0)
 #print('PSSM#', len(pssm),pssm)
 return(pssm)