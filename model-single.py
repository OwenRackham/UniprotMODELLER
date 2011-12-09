import sys
import os
sys.path.append('/usr/lib/python2.5/site-packages/')
from modeller import *
from modeller.automodel import *
sys.stdout = open(os.devnull,"w")
sys.stderr = open(os.devnull,"w")
os.chdir('../temp/spiricoil')


env = environ()
a = automodel(env, alnfile=sys.argv[1],
              knowns=sys.argv[2], sequence=sys.argv[3],
              assess_methods=(assess.DOPE, assess.GA341))
a.starting_model = 1
a.ending_model = 1
a.make()
os.chdir('../../cgi-bin')
