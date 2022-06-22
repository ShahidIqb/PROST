#!/usr/bin/env python
from misc_yang import *

rnam3_std0 = ["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU",
				"MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR", "UNK"]
rnam1_std0 = "ACDEFGHIKLMNPQRSTVWYX"
rnam3_std = rnam3_std0[:20]
rnam1_std = rnam1_std0[:20]
ASA_std = (115, 135, 150, 190, 210, 75, 195, 175, 200, 170,
			185, 160, 145, 180, 225, 115, 140, 155, 255, 230)
aa1_std = rnam1_std
aa3_std = rnam3_std
dict_rnam3_std = dict(zip(rnam3_std, range(len(rnam3_std))))
dict_rnam1_std = dict(zip(rnam1_std, range(len(rnam1_std))))
dict_rnam1_ASA = dict(zip(rnam1_std, ASA_std))
#
def aaDefine(rn, itype=0, DEBUG=0):
	if itype == 0: itype = len(rn)
	if itype == 3:
		ir = dict_rnam3_std.get(rn, -1)
	elif itype == 1:
		if rn == 'Z': rn = 'E'
		elif rn == 'B': rn = 'D'
		ir = rnam1_std.find(rn)
	else: die("not known residue type: %d %s" % (itype, rn))
	if DEBUG and ir < 0:
		print >>sys.stderr, 'unknown res:', rn
	return ir
def ssDefine(ss1):
	if ss1 in ' CST': return 0
	elif ss1 in 'HGI': return 1
	elif ss1 in 'EBb': return 2
	print >>sys.stderr, 'unknown ssec1: ', ss1
	sys.exit()
def rdfasta(fn):
	seq1 = ""; name = ''
	for i,x in enumerate(open(fn)):
		if i==0 and x[0]=='>': name = x[1:].strip()
		elif x[0] == '>': break
		else: seq1 += x.strip()
	return name, seq1
def rdseqs(fn, bstrip=1):
	list_name, list_seq = [], []
	for i,x in enumerate(open(fn)):
		if x[0] == '>':
			list_name.append(x[1:].strip())
			list_seq.append('')
		else:
			if len(list_seq) == 0 and x.strip()!='':
				list_name.append('UNK'); list_seq.append('')
#				print >>stderr, 'skip wrong rec: %s' % x.strip()
#				continue
			if bstrip: x2 = x.strip()
			else: x2 = x
			list_seq[-1] += x2
	return list_name, list_seq
def reverse_seq(seq1):
	dict_rev_seq = dict(zip('ATGCU', 'TACGA'))
	return ''.join([dict_rev_seq.get(a,'X') for a in seq1.upper()[::-1]])
def translate_dna(sequence):
	gencode = { 'ATA':'I',
	'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T',
	'ACT':'T', 'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S',
	'AGT':'S', 'AGA':'R', 'AGG':'R', 'CTA':'L', 'CTC':'L', 'CTG':'L',
	'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H',
	'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R',
	'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A',
	'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D', 'GAA':'E',
	'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 'TCA':'S',
	'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L',
	'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C',
	'TGT':'C', 'TGA':'*', 'TGG':'W', }
	sequence = sequence.upper().replace('U', 'T')
	if len(sequence) % 3 != 0:
		print >>stderr, 'seqlen is not 3s: %d' % len(sequence)
	proteinseq = ''
	for n in range(0,len(sequence),3):
		d1 = sequence[n:n+3]
		if gencode.has_key(d1) == True:
			proteinseq += gencode[d1]
		else: print >>stderr, ('unknown code: %s' % d1)
	return proteinseq
#
###
class SeqRec:
	def __init__(self):
		self.name = self.description = None
		self.seq = ''
	def __str__(self):
		return '%s| %s: %s ..' %(self.name, self.description, self.seq[:30])
#
class Parse_seq:
	def __init__(self, fn, bstrip=True, ipos=None, sep1=r'[>|:\s]+'):
		self.fp = open(fn)
		self.lastline = self.fp.next()
		self.bstrip = bstrip
		self.sep1 = sep1
		self.ipos = ipos
	def __iter__(self):
		return self
	def proc_seq(self, f1):
		if self.bstrip: return f1.rstrip()
		else: return f1
	def next(self):
		rec1 = SeqRec()
# minor change on 11/19/2015
# parse the head
		v1 = []
		if self.lastline is None: raise StopIteration
		elif self.lastline[0] == '>':
			rec1.description = info1 = self.lastline[1:].rstrip()
			rec1.name = ''
			if self.ipos is not None:
				ss = re.split(self.sep1, info1)
				if ss[0] == '': del ss[0]
				if len(ss) > self.ipos: rec1.name = ss[self.ipos]
			else:
				ss = info1.split()
				if len(ss) > 0: rec1.name = ss[0]
		else:
			v1 = [self.proc_seq(self.lastline)]
#
# read seq
		for x in self.fp:
			if x[0] == '>':
				self.lastline = x
				break
			else:
				v1.append( self.proc_seq(x) )
		else: self.lastline = None
		rec1.seq = ''.join(v1)

		return rec1
###
#
if __name__ == '__main__':
	for x in Parse_seq(argv[1]): print x
