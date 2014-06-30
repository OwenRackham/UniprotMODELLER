import sys
import os
sys.path.append('/usr/lib/python2.5/dist-packages/modeller/')
from modeller import *
from modeller.automodel import *
if __name__ == '__main__':
	#sys.stdout = open(os.devnull,"w")
	#sys.stderr = open(os.devnull,"w")
	path = '../data/%s' %(sys.argv[4])
	if os.path.exists(path):
		os.chdir(path)
	else:
		 os.makedirs(path)
		 os.chdir(path)
		 

	env = environ()
	a = automodel(env, alnfile=sys.argv[1],
	              knowns=sys.argv[2], sequence=sys.argv[3],
	              assess_methods=(assess.DOPE, assess.GA341))
	a.starting_model = 1
	a.ending_model = 1
	a.make()
	os.chdir('/home/rackham/workspace/UniprotMODELLER/')
