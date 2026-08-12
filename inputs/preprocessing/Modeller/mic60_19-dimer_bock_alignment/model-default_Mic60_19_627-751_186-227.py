from modeller import *
from modeller.automodel import *
log.verbose()
env = Environ() # create a new MODELLER environment to build this model in
# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']
a = AutoModel(env,
alnfile = 'Mic60_19.ali', # alignment filename
knowns = '7pv1_split',# codes of the templates
sequence = 'Mic60_19',
assess_methods=(assess.DOPE)) # code of the target
a.starting_model= 1 # index of the first model
a.ending_model = 20 # index of the last model
# (determines how many models to calculate)
a.make() # do the actual comparative modeling
