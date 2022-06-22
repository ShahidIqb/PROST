#!/usr/bin/env python
import sys, os, numpy, re
import numpy as np
from sys import *
#from math import *
""" Some general functions for PYTHON """

def open(fn, action='r'):
	""" open one file and return file pointer. Supported include:
	".gz" for gzip files
	"--" for stdin/stdout
	"""
	if isinstance(fn, file): return fn
	if 'w' in action:
		if fn in ['--', 'STDOUT']: return stdout
	else:
		if fn in ['--', 'STDIN']: return stdin
		elif fn.endswith('.gz'):
			import gzip
			return gzip.open(fn)
	return file(fn, action)
openfile = open
#
def rdlist_y(fn, ipos=0):
	""" read one file and return the first field as a list """
	list1 = []
	for x in open(fn):
		if x[0] == '#': continue
		ss = x.split()
		if len(ss) < 1: continue
		if ipos == 'ss': ss2 = ss
		elif isinstance(ipos, int): ss2 = ss[ipos]
		else: ss2 = [ss[k] for k in ipos]
		list1.append(ss2)
	return list1
#
def rddict_y(fn, ipos1=0, ipos2='ss'):
	""" read one file and return a dict of first field/line """
	dict1 = {}
	stype = ipos2
	if stype == 'line': stype = 1
	elif ipos2 == 'ss': stype = 2
	else:
		stype = 0
		try: ipos2 = int(ipos2)
		except:
			ipos2 = map(int, ipos2)
			stype = -1

	for x in open(fn):
		if x[0] == '#': continue
		ss = x.split()
		if len(ss) < 1: continue
		try:
			if ss[ipos1] in dict1:
				print >>stderr, 'duplicate key in 2 lines: ', ss[ipos1]
				print >>stderr, x, dict1[ss[ipos1]], '\n'
				continue
			if stype==0: v1 = ss[ipos2]
			elif stype==1: v1 = x.rstrip()
			elif stype==2: v1 = ss
			else: v1 = [ss[k] for k in ipos2]
		except: continue
		dict1[ss[ipos1]] = v1
	return dict1
#
def die(*info):
	"""  print info and exit """
	print >>stderr, list2str(info)
	exit(1)
#
def rdelete(line, exts):
	""" return string by deleting 'ext' in the right """
	if isinstance(exts, str): exts = [exts]
	for ext in exts:
		if line.endswith(ext): return line[:-len(ext)]
	return line
#
def list2str(l1, fmt=None, sep=' '):
	if fmt == None:
		return sep.join([str(x) for x in l1])
	else:
		return sep.join([fmt % x for x in l1])
#
def getfile(list_dir, list_bn, DEBUG=True):
	""" find files with existed combinations """
	for x in list_dir:
		for y in list_bn:
			fn = os.path.join(x, y)
			if os.path.isfile(fn): return fn
	if DEBUG: print >>stderr, 'no file found', list_dir, list_bn
	return ""
def fold_str(s1, n0=80):
	return '\n'.join([s1[k:k+n0] for k in range(0, len(s1), n0)])
#
#
#
def readopt(args, opt1, nout=1):
	''' simple way to read opt '''
	print >>stderr, "obsolete readopt (using find)"
	if opt1 not in args: return
	i0 = args.index(opt1)
	if len(args) <= nout + i0:
		die ('wrong opt: "%s" with %d paramters' % (opt1, nout))
		return
	elif nout == 1:
		return args[i0+1]
	else:
		return args[i0+1:i0+nout+1]
