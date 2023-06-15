
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
import math
import os
import sys

runType = sys.argv[1] # Specify test or prod
runID = sys.argv[2]   # Specify the number of runs
run_output_dir = 'run_' + str(runID)

if runType == "test":
    num_frames = 1000
elif runType == "prod":
    num_frames = 20000

max_shuffle_core = 5
max_rotation_core = 0.5238      # 30 degree
max_shuffle_set2 = 75
rex_max_temp = 1.25

#xlinkdata files
xl_BDP_PIR_human = '/home/muskaan/Documents/modeling_micos_complex/Data/crosslinks/human/BDP_PIR_human.csv'
xl_BDP_PIR_mouse = '/home/muskaan/Documents/modeling_micos_complex/Data/crosslinks/human/BDP_PIR_mouse.csv'
xl_DSS_BS3 = '/home/muskaan/Documents/modeling_micos_complex/Data/crosslinks/human/DSS_BS3_human.csv'
xl_DHSO_DSSO = '/home/muskaan/Documents/modeling_micos_complex/Data/crosslinks/human/DHSO_DSSO_human.csv' # yu. bartolec, xlinkdb dsso xlinks
xl_DSSO_mouse = '/home/muskaan/Documents/modeling_micos_complex/Data/crosslinks/human/DSSO_XLinkDB_mouse.csv'


# Topology File
topology_file = "/home/muskaan/Documents/modeling_micos_complex/Data/topology_ideal_helix.txt"



### --------------------------------
### WRAPPERS ###
### --------------------------------

# wrapper for the Biochemical Binding Restraint
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


# wrapper for Trans membrane Restraint
class TransMembraneRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist, R, r, sigma, weight=1):
        particles = plist
        name = 'TransMembraneRestraint%1%'
        super(TransMembraneRestraintC, self).__init__(model, name=name, weight=weight)
        res_main = IMP.micos.TransMembraneRestraint(particles, R, r, sigma)
        self.rs.add_restraint(res_main)
        print("Trans membrane Restraint applied")

# wrapper for surface localization restraint
class SurfaceLocalizationRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist,r, sigma, allowed_dist, weight=1):
        particles = plist
        name = 'SurfaceLocalizationRestraint%1%'
        super(SurfaceLocalizationRestraintC, self).__init__(model, name=name, weight=weight)
        res_main = IMP.micos.SurfaceLocalizationRestraint(particles, r, sigma, allowed_dist)
        self.rs.add_restraint(res_main)


# wrapper for IMS localization restraint
class IMSLocalizationRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist,r, sigma, weight=1):
        particles = plist
        name = 'IMSLocalizationRestraint%1%'
        super(IMSLocalizationRestraintC, self).__init__(model, name=name, weight=weight)
        res_main = IMP.micos.IMSLocalizationRestraint(particles, r, sigma)
        self.rs.add_restraint(res_main)

# wrapper for matrix localization restraint
class MatrixLocalizationRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist,R, sigma, weight=1):
        particles = plist
        name = 'MatrixLocalizationRestraint%1%'
        super(MatrixLocalizationRestraintC, self).__init__(model, name=name, weight=weight)
        res_main = IMP.micos.MatrixLocalizationRestraint(particles, R, sigma)
        self.rs.add_restraint(res_main)

# wrapper for z axial Restraint
class ZAxialRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist, cj, om, sigma, weight=1):
        particles = plist
        name = 'ZAxialRestraint%1%'
        super(ZAxialRestraintC, self).__init__(model, name=name, weight=weight)
        res_main = IMP.micos.ZAxialRestraint(particles, cj, om, sigma)
        self.rs.add_restraint(res_main)
        print("zaxial Restraint applied")


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
root_hier, dof = bs.execute_macro(max_rb_trans= 0.73,
                                  max_rb_rot= 0.23,
                                  max_bead_trans= 3.29,
                                  max_srb_trans= 0.1,
                                  max_srb_rot=0.01)

# molecules = t.get_components()
# print(molecules)
# exit()

# Uncomment the following lines to get test.rmf file to visualise the system representation

# Uncomment this line for verbose output of the representation
# IMP.atom.show_with_representations(root_hier)
# # output to RMF
# fname = 'test.rmf'
# rh = RMF.create_rmf_file(fname)
# IMP.rmf.add_hierarchy(rh, root_hier)
# IMP.rmf.save_frame(rh)


# exit()
################################################################################
########################## Fixing Particles ####################################
################################################################################
# # First select and gather all particles to fix.
# mem_surface_particles = [(634, 637,'MIC60', None, None),(643, 644,'MIC60', None, None),(693,693,'MIC60', None, None)]
# # fixed_particles = IMP.pmi.tools.select_by_tuple_2(root_hier,mem_surface_particles,10)
# fixed_particles = []
# fixed_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC60',copy_indexes = [0,1],residue_indexes = range(627,752)).get_selected_particles()
# fixed_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC19',residue_indexes = range(186,227)).get_selected_particles()
#
# print(fixed_particles)
#
# # # Fix the Corresponding Rigid movers and Super Rigid Body movers using dof
# # # The flexible beads will still be flexible (fixed_beads is an empty list)!
# fixed_beads,fixed_rbs=dof.disable_movers(fixed_particles, [IMP.core.RigidBodyMover,
#                                          IMP.pmi.TransformMover])
# # # print('######################################')
# print(fixed_rbs)
# # print(fixed_beads)
# # exit()

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
#
# # # -------------------------------
# # MEMBRANE RESTRAINTS
# these restraints are applied to localize the beads w.r.t. membrane topology
# there are 4 restaints: Trans membrane for TM regions, IMS localization to localize beads within the inner radius,
# Matrix localization to keep beads farther away from outer radius and surface localization to keep the beads close to inner cylinder
# by defining allowed_dist and keeping the beads between the allowed_dist and the inner radius
# add center beads for membrane surface so that it can score based on those beads only
# a cylinder can be used (using BILD files) to visualize the restraints
#
# select_by_tuple_2(start,stop,molname,copynum,statenum);  use 'None' for them which will get all copies and states
# # # -------------------------------

ims_particles = []
tm_particles = []
mem_surface_particles = []
matrix_particles = []
above_cj_particles = []

ims_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC60',copy_indexes = [0,1],residue_indexes = range(410,752)).get_selected_particles()
ims_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC19',residue_indexes = range(1,227)).get_selected_particles()
ims_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC60',copy_indexes = [2,3],residue_indexes = range(172,582)).get_selected_particles()
ims_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC13',residue_indexes = range(24,118)).get_selected_particles()
ims_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC10',residue_indexes = range(1,12)).get_selected_particles()
ims_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC10',residue_indexes = range(61,78)).get_selected_particles()

tm_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC10',residue_indexes = range(13,36)).get_selected_particles()
tm_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC10',residue_indexes = range(40,60)).get_selected_particles()
tm_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC60',residue_indexes = range(149,171)).get_selected_particles()
tm_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC13',residue_indexes = range(8,23)).get_selected_particles()

matrix_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC13',residue_indexes = range(1,7)).get_selected_particles()
matrix_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC10',residue_indexes = range(37,39)).get_selected_particles()
matrix_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC60',residue_indexes = range(1,148)).get_selected_particles()

mem_surface_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC60',residue_indexes = range(627,758)).get_selected_particles()
mem_surface_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC19',residue_indexes = range(186,227)).get_selected_particles()

above_cj_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC19',residue_indexes = range(1,50)).get_selected_particles()


output_objects = []
R = 130 #radius of the outer cylinder
r = 90 #radius of the inner cylinder
sigma = 0.002 #weight
allowed_dist = 5 #for surface localization, this is the distance from the center. the beads should stay within this and inner radius
cj = 0
om = -50

mir = TransMembraneRestraintC(mdl,tm_particles, R,r,sigma)
output_objects.append(mir)
mir.add_to_model()
print("Trans membrane Restraint applied")


slr = SurfaceLocalizationRestraintC(mdl,mem_surface_particles, r, sigma, allowed_dist)
output_objects.append(slr)
slr.add_to_model()
print("surface localization restraint applied")


ilr = IMSLocalizationRestraintC(mdl,ims_particles,r,sigma)
output_objects.append(ilr)
ilr.add_to_model()
print("IMS localization restraint applied")


mlr = MatrixLocalizationRestraintC(mdl,matrix_particles,R,sigma)
output_objects.append(mlr)
mlr.add_to_model()
print("matrix localization restraint applied")

zar = ZAxialRestraintC(mdl,above_cj_particles,cj,om,sigma)
output_objects.append(zar)
zar.add_to_model()
print("zaxial restraint applied")

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
#
xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
xldbkc.set_standard_keys()

xldb_BDP_PIR_human = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_BDP_PIR_human.create_set_from_file(file_name=xl_BDP_PIR_human,
                                 converter=xldbkc)
xlr_BDP_PIR_human = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                database=xldb_BDP_PIR_human, # The crosslink database.
                length=33,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.02,           # This adds a linear term to the scoring function
                label="BDP_PIR_human",                        #   to bias crosslinks towards each other
                weight=10)       # Scaling factor for the restraint score.

output_objects.append(xlr_BDP_PIR_human)
xlr_BDP_PIR_human.add_to_model()

xldb_BDP_PIR_mouse = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_BDP_PIR_mouse.create_set_from_file(file_name=xl_BDP_PIR_mouse,
                                 converter=xldbkc)
xlr_BDP_PIR_mouse = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                database=xldb_BDP_PIR_mouse, # The crosslink database.
                length=33,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.02,           # This adds a linear term to the scoring function
                label="BDP_PIR_mouse",                        #   to bias crosslinks towards each other
                weight=10)       # Scaling factor for the restraint score.

output_objects.append(xlr_BDP_PIR_mouse)
xlr_BDP_PIR_mouse.add_to_model()

xldb_DSS_BS3 = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_DSS_BS3.create_set_from_file(file_name=xl_DSS_BS3,
                                 converter=xldbkc)
xlr_DSS_BS3 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                database=xldb_DSS_BS3, # The crosslink database.
                length=25,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.02,           # This adds a linear term to the scoring function
                label="DSS_BS3",                        #   to bias crosslinks towards each other
                weight=10)       # Scaling factor for the restraint score.

output_objects.append(xlr_DSS_BS3)
xlr_DSS_BS3.add_to_model()

xldb_DHSO_DSSO = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_DHSO_DSSO.create_set_from_file(file_name=xl_DHSO_DSSO,
                                 converter=xldbkc)
xlr_DHSO_DSSO = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                database=xldb_DHSO_DSSO, # The crosslink database.
                length=21,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.02,           # This adds a linear term to the scoring function
                label="DHSO_DSSO",                        #   to bias crosslinks towards each other
                weight=10)       # Scaling factor for the restraint score.

output_objects.append(xlr_DHSO_DSSO)
xlr_DHSO_DSSO.add_to_model()

xldb_DSSO_mouse = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_DSSO_mouse.create_set_from_file(file_name=xl_DSSO_mouse,
                                 converter=xldbkc)
xlr_DSSO_mouse = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
                root_hier=root_hier,    # Must pass the root hierarchy to the system
                database=xldb_DSSO_mouse, # The crosslink database.
                length=21,              # The crosslinker plus side chain length
                resolution=1,           # The resolution at which to evaluate the crosslink
                slope=0.02,           # This adds a linear term to the scoring function
                label="DSSO_mouse",                        #   to bias crosslinks towards each other
                weight=10)       # Scaling factor for the restraint score.

output_objects.append(xlr_DSSO_mouse)
xlr_DSSO_mouse.add_to_model()

# -----------------------------
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
evr.add_to_model()

print("Excluded volume restraint applied")




# #####################################################
# ###################### SAMPLING #####################
# #####################################################
# With our representation and scoring functions determined, we can now sample
# the configurations of our model with respect to the information.
# print("The type of run is: " + str(runType))
# print("Number of sampling frames: " + str(num_frames))


############ SHUFFLING ##################
## 3 bounding boxes for 3 regions - membrane, IMS and Matrix to shuffle the proteins in that region only

ims_bb = ((0,-r,0),(r,0,100))
tm_bb = ((90,90,0),(130,130,100))

# First shuffle all particles to randomize the starting point of the
# system. For larger systems, you may want to increase max_translation

IMP.pmi.tools.shuffle_configuration(ims_particles,
                                    max_translation=50,
                                    bounding_box = ims_bb,
                                    avoidcollision_rb = False)


IMP.pmi.tools.shuffle_configuration(tm_particles,
                                    max_translation=50,
                                    bounding_box = tm_bb,
                                    avoidcollision_rb = False)

IMP.pmi.tools.shuffle_configuration(matrix_particles,
                                    max_translation=50)

# Shuffling randomizes the bead positions. It's good to
# allow these to optimize first to relax large connectivity
# restraint scores.  100-500 steps is generally sufficient.
dof.optimize_flexible_beads(1000)


IMP.pmi.dof.DegreesOfFreedom.enable_all_movers(dof)

print("Replica Exchange Maximum Temperature : " + str(rex_max_temp))

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
        root_hier=root_hier,                    # pass the root hierarchy
        crosslink_restraints=[xlr_BDP_PIR_human, xlr_BDP_PIR_mouse, xlr_DSS_BS3, xlr_DHSO_DSSO, xlr_DSSO_mouse],
        # This allows viewing the crosslinks in Chimera. Also, there is not inter-protein ADH crosslink available. Hence it is not mentioned in this list
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
