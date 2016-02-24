from modeller import *
from modeller.automodel import *
#from modeller import soap_protein_od

env = environ()
a = automodel(env, alnfile='OR414-2ctf.ali',
              knowns='2ctf', sequence='OR414',
              assess_methods=(assess.DOPE,
                              # soap_protein_od.Scorer(),
                              assess.GA341))
a.starting_model = 1
a.ending_model = 5
a.make()
