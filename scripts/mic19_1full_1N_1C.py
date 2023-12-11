
##########################################################################################
######################### IMP Modeling Script for MICOS Complex ##########################
##########################################################################################


# Imports
from __future__ import print_function
import IMP
import RMF
import IMP.rmf
import IMP.pmi
import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.topology
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.dof
import IMP.atom
import IMP.micos
import IMP.container
import math
import os
import sys
import ihm


# ### --------------------------------
### WRAPPERS ###
### --------------------------------

# Wrapper for the Biochemical Binding Restraint
class MinimumPairDistanceBindingRestraint(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist1, plist2, x0=0, kappa=1, label=None, weight=1):
        name = 'MinimumPairDistanceBindingRestraint%1%'
        super(MinimumPairDistanceBindingRestraint, self).__init__(model, name=name, label=label, weight=weight)
        l1 = IMP.container.ListSingletonContainer(mdl)
        l1.add(plist1)
        l2 = IMP.container.ListSingletonContainer(mdl)
        l2.add(plist2)
        bipartite_container = IMP.container.AllBipartitePairContainer(l1, l2)
        score = IMP.core.HarmonicSphereDistancePairScore(x0, kappa)
        res_main = IMP.container.MinimumPairRestraint(score, bipartite_container, 1)
        self.rs.add_restraint(res_main)
        print("RESTRAINT: Added MPDBR on particles", len(plist1), ":", len(plist2), "at x0:kappa", str(x0), ":",
              str(kappa))

# Note: sigma is mem_wt
# Wrapper for Transmembrane Restraint
class TransMembraneRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist, R, r, sigma, weight=1):
        particles = plist
        name = 'TransMembraneRestraint%1%'
        super(TransMembraneRestraintC, self).__init__(model, name=name, weight=weight)
        res_main = IMP.micos.TransMembraneRestraint(particles, R, r, sigma)
        self.rs.add_restraint(res_main)


# Wrapper for surface localization restraint
class SurfaceLocalizationRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist,r, sigma, allowed_dist, weight=1):
        particles = plist
        name = 'SurfaceLocalizationRestraint%1%'
        super(SurfaceLocalizationRestraintC, self).__init__(model, name=name, weight=weight)
        res_main = IMP.micos.SurfaceLocalizationRestraint(particles, r, sigma, allowed_dist)
        self.rs.add_restraint(res_main)


# Wrapper for IMS localization restraint
class IMSLocalizationRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist,r, sigma, weight=1):
        particles = plist
        name = 'IMSLocalizationRestraint%1%'
        super(IMSLocalizationRestraintC, self).__init__(model, name=name, weight=weight)
        res_main = IMP.micos.IMSLocalizationRestraint(particles, r, sigma)
        self.rs.add_restraint(res_main)

# Wrapper for matrix localization restraint
class MatrixLocalizationRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist,R, sigma, weight=1):
        particles = plist
        name = 'MatrixLocalizationRestraint%1%'
        super(MatrixLocalizationRestraintC, self).__init__(model, name=name, weight=weight)
        res_main = IMP.micos.MatrixLocalizationRestraint(particles, R, sigma)
        self.rs.add_restraint(res_main)

# Wrapper for zaxial restraint on protein level
class ZAxialRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist, ub, lb, sigma, method,label=None, weight=1):
        particles = plist
        name = 'ZAxialRestraint%1%'
        super(ZAxialRestraintC, self).__init__(model, name=name, weight=weight, label=label)
        res_main = IMP.micos.ZAxialRestraint(particles, ub, lb, sigma,method)
        self.rs.add_restraint(res_main)


runType = sys.argv[1] # Specify test or prod
runID = sys.argv[2]   # Specify the number of runs
run_output_dir = 'run_' + runID
data_direc = sys.argv[3]

if runType == "test":
    num_frames = 1500
elif runType == "prod":
    num_frames = 30000



#Xlinkdata files
xl_BDP_PIR_human = f'{data_direc}/crosslinks/human/sampling_BDP_PIR_human'
xl_BDP_PIR_mouse = f'{data_direc}/crosslinks/human/sampling_BDP_PIR_mouse'
xl_DSS_BS3_human = f'{data_direc}/crosslinks/human/sampling_DSS_BS3_human'
xl_DHSO_DSSO_human = f'{data_direc}/crosslinks/human/sampling_DHSO_DSSO_human' # yu. bartolec, xlinkdb dsso xlinks
xl_DSSO_mouse = f'{data_direc}/crosslinks/human/sampling_DSSO_mouse'


# Topology File
topology_file = f'{data_direc}/topology_mic19_1full_1N_1C_finer_rep.txt'

# Weights
mem_wt = 0.04
zar_wt = 0.2
coIP_wt = 4
BNPAGE_wt = 2
WB_wt = 1
AF_wt = 5


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Here is where the work begins
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# All IMP systems start out with a Model
mdl = IMP.Model()

# Read the topology file for a given state
t = IMP.pmi.topology.TopologyReader(topology_file)


# Create a BuildSystem macro to and add a state from a topology file
bs = IMP.pmi.macros.BuildSystem(mdl)

bs.add_state(t)

# executing the macro will return the root hierarchy and degrees of freedom (dof) objects
root_hier, dof = bs.execute_macro(max_rb_trans= 0.1,
                                  max_rb_rot= 0.1,
                                  max_bead_trans= 2.43,
                                  max_srb_trans= 0.1,
                                  max_srb_rot= 0.02)

rex_max_temp = 1.25

# Uncomment the following lines to get test.rmf file to visualise the system representation

# Uncomment this line for verbose output of the representation
# IMP.atom.show_with_representations(root_hier)
# # output to RMF
# fname = 'test.rmf'
# rh = RMF.create_rmf_file(fname)
# IMP.rmf.add_hierarchy(rh, root_hier)
# IMP.rmf.save_frame(rh)
#
# exit()
#####################################################
##################### RESTRAINTS ####################
#####################################################

# Restraints define functions that score the model based on
# input information.
#
# Restraint objects are first created in the definition.
# To be evaluated, the restraint object must be add_to_model().
#
# In some cases, sampled parameters for restraints must be added to the DOF
# object

# The output_objects list is used to collect all restraints
# where we want to log the output in the STAT file.
# Each restraint should be appended to this list.

output_objects = []

# # # -------------------------------
# # MEMBRANE RESTRAINTS
# these restraints are applied to localize the beads w.r.t. membrane topology
# there are 4 restaints: Trans membrane for TM regions, IMS localization to localize beads within the inner radius,
# Matrix localization to keep beads farther away from outer radius and surface localization to keep the beads close to inner cylinder
# by defining allowed_dist and keeping the beads between the allowed_dist and the inner radius
# add center beads for membrane surface so that it can score based on those beads only
# a cylinder can be used (using BILD files) to visualize the restraints
# # # -------------------------------

## These selections are for the membrane restraints.
mic10_selections = []
mic10_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC10',residue_indexes = range(1,13)).get_selected_particles()) # ims
mic10_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC10',residue_indexes = range(13,37)).get_selected_particles()) # tm
mic10_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC10',residue_indexes = range(37,40)).get_selected_particles()) # matrix
mic10_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC10',residue_indexes = range(40,61)).get_selected_particles()) # tm
mic10_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC10',residue_indexes = range(61,79)).get_selected_particles()) # ims

mic60_selections = []
mic60_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC60',copy_indexes = [0,1],residue_indexes = range(627,759)).get_selected_particles()) #ims
mic60_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC60',copy_indexes = [0,1],residue_indexes = range(583,627)).get_selected_particles()) # above z
mic60_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC60',copy_indexes = [0,1,2,3],residue_indexes = range(410,583)).get_selected_particles()) # on the center
mic60_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC60',copy_indexes = [0,1],residue_indexes = range(634,636)).get_selected_particles()) # slr
mic60_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC60',copy_indexes = [0,1],residue_indexes = range(643,645)).get_selected_particles()) # slr
# mic60_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC60',copy_indexes = [0,1],residue_indexes = range(692,693)).get_selected_particles()) # slr

mic13_selections = []
mic13_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC13',residue_indexes = range(1,8)).get_selected_particles()) # matrix
mic13_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC13',residue_indexes = range(8,24)).get_selected_particles()) # tm
mic13_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC13',residue_indexes = range(24,119)).get_selected_particles()) # ims

mic19_selections = []
mic19_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC19',copy_indexes = [0,1],residue_indexes = range(1,186)).get_selected_particles()) # ims
mic19_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC19',copy_indexes = [0,1],residue_indexes = range(1,15)).get_selected_particles()) # above z
mic19_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC19',copy_indexes = [0,2],residue_indexes = range(186,228)).get_selected_particles()) # above z
mic19_selections.append(IMP.atom.Selection(hierarchy = root_hier,molecule='MIC19',copy_indexes = [0,1],residue_indexes = range(59,175)).get_selected_particles()) # on the center

R = 130 #radius of the outer cylinder
r = 90 #radius of the inner cylinder
allowed_dist = 10 #for surface localization, this is the distance from the center. the beads should stay within this and the radius of the inner cylinder
lb = 0
ub = 50
max_in_cj = 100
min_up_cj = -70

# TransMembraneRestraint
for particle_set in (mic10_selections[1], mic10_selections[3], mic13_selections[1]):
    tmr = TransMembraneRestraintC(mdl,particle_set,R,r,mem_wt)
    output_objects.append(tmr)
    tmr.add_to_model()

# SurfaceLocalizationRestraint
for particle_set in (mic60_selections[3],mic60_selections[4]):
    slr = SurfaceLocalizationRestraintC(mdl,particle_set,r,mem_wt,allowed_dist)
    output_objects.append(slr)
    slr.add_to_model()

#IMSLocalizationRestraint
for particle_set in (mic10_selections[0],mic10_selections[4],mic13_selections[2],mic60_selections[2], mic60_selections[1],mic19_selections[2],mic60_selections[0],mic19_selections[0]):
    ilr = IMSLocalizationRestraintC(mdl,particle_set,r,mem_wt)
    output_objects.append(ilr)
    ilr.add_to_model()

#MatrixLocalizationRestraint
for particle_set in (mic10_selections[2],mic13_selections[0]):
    mlr = MatrixLocalizationRestraintC(mdl,particle_set,R,mem_wt)
    output_objects.append(mlr)
    mlr.add_to_model()

#ZAxialRestraint for interaction with SAM/TOB complex
## for Mic60 ateast one bead above cj
for particle_set in (mic60_selections[1], mic60_selections[0]):
    zar_sam50 = ZAxialRestraintC(mdl,particle_set,min_up_cj,lb,zar_wt,'domain',label='zar_sam50')
    output_objects.append(zar_sam50)
    zar_sam50.add_to_model()

## for Mic19 N-term, larger bounds
zar_sam50 = ZAxialRestraintC(mdl,mic19_selections[1],min_up_cj,ub,zar_wt,'None',label='zar_sam50_mic19')
output_objects.append(zar_sam50)
zar_sam50.add_to_model()

#ZAxialRestraint at the rim of the cylinder
zar_cc = ZAxialRestraintC(mdl,mic60_selections[2],lb,ub,zar_wt,'None',label = 'zar_cc')
output_objects.append(zar_cc)
zar_cc.add_to_model()

#ZAxialRestraint in the cylinder
for particle_set in (mic10_selections[0],mic10_selections[1],mic10_selections[2],mic10_selections[3],mic10_selections[4],mic13_selections[0],mic13_selections[1],mic13_selections[2]):
    zar_below = ZAxialRestraintC(mdl,particle_set,lb,max_in_cj,zar_wt,'average',label='zar_below')
    output_objects.append(zar_below)
    zar_below.add_to_model()


# To specify the particles in the specific bounding box
# ims_particles = []
# mem_surface_particles = []
tm_particles = []
matrix_particles = []
ims_mem_surface_particles = []

ims_mem_surface_particles.extend([mic10_selections[0], mic10_selections[4], mic13_selections[2], mic60_selections[2], mic60_selections[1], mic60_selections[0], mic19_selections[0], mic19_selections[2]])
tm_particles.extend([mic10_selections[1], mic10_selections[3], mic13_selections[1]])
matrix_particles.extend([mic10_selections[2],mic13_selections[0]])
# ims_particles.extend([mic10_selections[0], mic10_selections[4], mic13_selections[2], mic60_selections[2], mic60_selections[1], mic19_selections[0], mic19_selections[2]])
# mem_surface_particles.extend([mic60_selections[1],mic19_selections[2]])

# # -----------------------------
# %%%%% CONNECTIVITY RESTRAINT
#
# Restrains residues/particles that are connected in sequence
# This should be used for any system without an atomic force field (e.g. CHARMM)
# We apply the restraint to each molecule

for m in root_hier.get_children()[0].get_children():
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(m)
    cr.add_to_model()
    output_objects.append(cr)

print("Connectivity restraint applied")


# -----------------------------
# %%%%% EXCLUDED VOLUME RESTRAINT
#
# Keeps particles from occupying the same area in space.
# Here, we pass a list of all molecule chains to included_objects to apply this to every residue.
# We could also have passed root_hier to obtain the same behavior.
#
# resolution=1000 applies this expensive restraint to the lowest resolution for each particle.
evr = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
                                            included_objects=[root_hier],
                                            resolution=1000)
output_objects.append(evr)
print("Excluded volume restraint applied")



# #####################################################
# ###################### SAMPLING #####################
# #####################################################
# With our representation and scoring functions determined, we can now sample
# the configurations of our model with respect to the information.

############ SHUFFLING ##################
## 3 bounding boxes for 3 regions - membrane, IMS and Matrix to shuffle the proteins in that region only
ims_bb = ((0,-r,0),(r,r,100)) # this is for half cylinder
tm_bb = ((r,-r,0),(R,R,100))

# First shuffle all particles to randomize the starting point of the
# system. For larger systems, you may want to increase max_translation

IMP.pmi.tools.shuffle_configuration(ims_mem_surface_particles,
                                    max_translation=50,
                                    bounding_box = ims_bb,
                                    avoidcollision_rb = False)

# IMP.pmi.tools.shuffle_configuration(mem_surface_particles,
#                                     max_translation=50,
#                                     bounding_box = ims_bb,
#                                     avoidcollision_rb = False)

IMP.pmi.tools.shuffle_configuration(tm_particles,
                                    max_translation=50,
                                    bounding_box = tm_bb,
                                    avoidcollision_rb = False)
#
IMP.pmi.tools.shuffle_configuration(matrix_particles,
                                    max_translation=50)
# Shuffling randomizes the bead positions. It's good to
# allow these to optimize first to relax large connectivity
# restraint scores.  100-500 steps is generally sufficient.
dof.optimize_flexible_beads(1000)
IMP.pmi.dof.DegreesOfFreedom.enable_all_movers(dof)

print("Replica Exchange Maximum Temperature : " + str(rex_max_temp))


# -----------------------------
# %%%%% MINIMUM PAIR RESTRAINT
# co-IP > BN-PAGE > WB as 1, 2, 4, respectively
mpr1 = MinimumPairDistanceBindingRestraint(mdl,IMP.pmi.tools.select_by_tuple_2(root_hier,(81,85,"MIC13",None,None),10),IMP.pmi.tools.select_by_tuple_2(root_hier,(1,78,"MIC10",None,None),10),0,1,"MIC13_RDWN_MIC10_mpr",coIP_wt)
mpr2 = MinimumPairDistanceBindingRestraint(mdl,IMP.pmi.tools.select_by_tuple_2(root_hier,(24,28,"MIC10",None,None),10),IMP.pmi.tools.select_by_tuple_2(root_hier,(1,78,"MIC10",None,None),10),0,1,"MIC10_MIC10_mpr1",BNPAGE_wt)
mpr3 = MinimumPairDistanceBindingRestraint(mdl,IMP.pmi.tools.select_by_tuple_2(root_hier,(46,52,"MIC10",None,None),10),IMP.pmi.tools.select_by_tuple_2(root_hier,(1,78,"MIC10",None,None),10),0,1,"MIC10_MIC10_mpr2",BNPAGE_wt)
mpr4 = MinimumPairDistanceBindingRestraint(mdl,IMP.pmi.tools.select_by_tuple_2(root_hier,(2,26,"MIC13",None,None),10),IMP.pmi.tools.select_by_tuple_2(root_hier,(410,758,"MIC60",None,None),10),0,1,"MIC13_MIC60_mpr1",coIP_wt)


output_objects.append(mpr1)
output_objects.append(mpr2)
output_objects.append(mpr3)
output_objects.append(mpr4)

mpr1.add_to_model()
mpr2.add_to_model()
mpr3.add_to_model()
mpr4.add_to_model()

## AF-multimer -------------------------------
# mic10-13, mic10-60, mic13-60

AF1 = MinimumPairDistanceBindingRestraint(mdl,IMP.pmi.tools.select_by_tuple_2(root_hier,(18,22,"MIC13"),10),IMP.pmi.tools.select_by_tuple_2(root_hier,(23,27,"MIC10"),10),0,1,"MIC13_TM_MIC10_af1",AF_wt)
AF2 = MinimumPairDistanceBindingRestraint(mdl,IMP.pmi.tools.select_by_tuple_2(root_hier,(88,92,"MIC13"),10),IMP.pmi.tools.select_by_tuple_2(root_hier,(657,668,"MIC60"),10),0,1,"MIC13_MIC60_af2",AF_wt)
AF3 = MinimumPairDistanceBindingRestraint(mdl,IMP.pmi.tools.select_by_tuple_2(root_hier,(3,3,"MIC10"),10),IMP.pmi.tools.select_by_tuple_2(root_hier,(639,642,"MIC60"),10),0,1,"MIC10_MIC60_af3",AF_wt)
AF4 = MinimumPairDistanceBindingRestraint(mdl,IMP.pmi.tools.select_by_tuple_2(root_hier,(2,3,"MIC10"),10),IMP.pmi.tools.select_by_tuple_2(root_hier,(720,724,"MIC60"),10),0,1,"MIC10_MIC60_af4",AF_wt)

output_objects.append(AF1)
output_objects.append(AF2)
output_objects.append(AF3)
output_objects.append(AF4)

AF1.add_to_model()
AF2.add_to_model()
AF3.add_to_model()
AF4.add_to_model()

# -------------------------
# %%%%% CROSSLINKING RESTRAINT
#
# Restrains two particles via a distance restraint based on
# an observed crosslink.
#
# First, create the crosslinking database from the input file
# The "standard keys" correspond to a crosslink csv file of the form:
#
# Protein1,Residue1,Protein2,Residue2
# A,18,G,24
# A,18,G,146
# A,50,G,146
# A,50,G,171
# A,50,G,189
#
# This restraint allows for ambiguity in the crosslinked residues,
# a confidence metric for each crosslink and multiple states.
# See the PMI documentation or the MMB book chapter for a
# full discussion of implementing crosslinking restraints.

# This first step is used to translate the crosslinking data file.
# The KeywordsConverter maps a column label from the xl data file
# to the value that PMI understands.
# Here, we just use the standard keys.
# One can define custom keywords using the syntax below.
# For example if the Protein1 column header is "prot_1"
# xldbkc["Protein1"]="prot_1"

# The CrossLinkDataBase translates and stores the crosslink information
# from the file "xl_data" using the KeywordsConverter.

xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkc.set_standard_keys()

xldb_BDP_PIR_human = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_BDP_PIR_human.create_set_from_file(file_name=xl_BDP_PIR_human,
                                 converter=xldbkc)
xlr_BDP_PIR_human = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                database=xldb_BDP_PIR_human, # The crosslink database.
                length=42,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.0001,           # This adds a linear term to the scoring function
                label="BDP_PIR_human",                        #   to bias crosslinks towards each other
                weight=10,                              # Scaling factor for the restraint score.
                linker=ihm.ChemDescriptor("bruce"))

output_objects.append(xlr_BDP_PIR_human)


xldb_BDP_PIR_mouse = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_BDP_PIR_mouse.create_set_from_file(file_name=xl_BDP_PIR_mouse,
                                 converter=xldbkc)
xlr_BDP_PIR_mouse = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                database=xldb_BDP_PIR_mouse, # The crosslink database.
                length=42,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.0001,           # This adds a linear term to the scoring function
                label="BDP_PIR_mouse",                        #   to bias crosslinks towards each other
                weight=8,                                     # Scaling factor for the restraint score.
                linker=ihm.ChemDescriptor("bruce"))

output_objects.append(xlr_BDP_PIR_mouse)


xldb_DSS_BS3 = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_DSS_BS3.create_set_from_file(file_name=xl_DSS_BS3_human,
                                 converter=xldbkc)
xlr_DSS_BS3 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                database=xldb_DSS_BS3, # The crosslink database.
                length=25,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.0001,           # This adds a linear term to the scoring function
                label="DSS_BS3",                        #   to bias crosslinks towards each other
                weight=10,                              # Scaling factor for the restraint score.
                linker=ihm.ChemDescriptor("bruce"))

output_objects.append(xlr_DSS_BS3)



xldb_DHSO_DSSO = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_DHSO_DSSO.create_set_from_file(file_name=xl_DHSO_DSSO_human,
                                 converter=xldbkc)
xlr_DHSO_DSSO = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                database=xldb_DHSO_DSSO, # The crosslink database.
                length=21,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.0001,           # This adds a linear term to the scoring function
                label="DHSO_DSSO",                        #   to bias crosslinks towards each other
                weight=10,                              # Scaling factor for the restraint score.
                linker=ihm.ChemDescriptor("bruce"))

output_objects.append(xlr_DHSO_DSSO)


xldb_DSSO_mouse = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_DSSO_mouse.create_set_from_file(file_name=xl_DSSO_mouse,
                                 converter=xldbkc)
xlr_DSSO_mouse = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                database=xldb_DSSO_mouse, # The crosslink database.
                length=21,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.0001,           # This adds a linear term to the scoring function
                label="DSSO_mouse",                        #   to bias crosslinks towards each other
                weight=8,                       # Scaling factor for the restraint score.
                linker=ihm.ChemDescriptor("bruce"))
output_objects.append(xlr_DSSO_mouse)

xlr_BDP_PIR_human.add_to_model()
xlr_BDP_PIR_mouse.add_to_model()
xlr_DSS_BS3.add_to_model()
xlr_DHSO_DSSO.add_to_model()
xlr_DSSO_mouse.add_to_model()

evr.add_to_model()
# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange(mdl,
        root_hier=root_hier,                    # pass the root hierarchy
        monte_carlo_temperature = 1.0,
        replica_exchange_minimum_temperature = 1.0,
        replica_exchange_maximum_temperature = rex_max_temp,
	    monte_carlo_sample_objects=dof.get_movers(),  # pass all objects to be moved ( almost always dof.get_movers() )
        global_output_directory=run_output_dir,      # The output directory for this sampling run.
        output_objects=output_objects,          # Items in output_objects write information to the stat file.
        monte_carlo_steps=10,                   # Number of MC steps between writing frames
        number_of_best_scoring_models=0,        # set >0 to store best PDB files (but this is slow)
        number_of_frames=num_frames)            # Total number of frames to run / write to the RMF file.
        #test_mode=test_mode)                    # (Ignore this) Run in test mode (don't write anything)

# Ok, now we finally do the sampling!
rex.execute_macro()
#
# Outputs are then analyzed in a separate analysis script.
