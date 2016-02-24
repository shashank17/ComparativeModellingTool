from modeller import *

env = environ()
aln = alignment(env)
mdl = model(env, file='4rrf', model_segment=('FIRST:A', 'LAST:A'))
aln.append_model(mdl, align_codes='4rrf', atom_files='4rrf.pdb')
aln.append(file='OR414.ali', align_codes='OR414')
aln.align2d()
aln.write(file='OR414-4rrf.ali', alignment_format='PIR')
aln.write(file='OR414-4rrf.pap', alignment_format='PAP')
