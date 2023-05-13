
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
import IMP.pmi.restraints.em
import IMP.pmi.dof
import IMP.atom
import IMP.micos
import math
# from fixed_shuffle import shuffle_configuration
import os
import sys

runType = sys.argv[1] # Specify test or prod
runID = sys.argv[2]   # Specify the number of runs
run_output_dir = 'run_' + str(runID)

if runType == "test":
    num_frames = 100
elif runType == "prod":
    num_frames = 20000

max_shuffle_core = 5
max_rotation_core = 0.5238      # 30 degree
max_shuffle_set2 = 75
rex_max_temp = 1.5


# xl_data = '../../../Data/inputs/xlinks/out_inter_xl.csv'
# Topology File
topology_file = "../micos/modeling/Data/topology.txt"



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


# wrapper for membrane inclusion restraint
class MembraneInclusionRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist, R, r, sigma, weight=1):
        particles = plist
        name = 'MembraneInclusionRestraint%1%'
        super(MembraneInclusionRestraintC, self).__init__(model, name=name, weight=weight)
        res_main = IMP.micos.MembraneInclusionRestraint(particles, R, r, sigma)
        self.rs.add_restraint(res_main)
        # print("Membrane inclusion Restraint applied")

# wrapper for surface localization restraint
class SurfaceLocalizationRestraintC(IMP.pmi.restraints.RestraintBase):

    def __init__(self, model, plist,r, sigma, max_limit, weight=1):
        particles = plist
        name = 'SurfaceLocalizationRestraint%1%'
        super(SurfaceLocalizationRestraintC, self).__init__(model, name=name, weight=weight)
        res_main = IMP.micos.SurfaceLocalizationRestraint(particles, r, sigma, max_limit)
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
root_hier, dof = bs.execute_macro(max_rb_trans= 1,
                                  max_rb_rot= 0.1,
                                  max_bead_trans= 3.2,
                                  max_srb_trans= 0.01,
                                  max_srb_rot=0.04)




molecules = t.get_components()


# Uncomment the following lines to get test.rmf file to visualise the system representation

# Uncomment this line for verbose output of the representation
IMP.atom.show_with_representations(root_hier)
# output to RMF
fname = 'test.rmf'
rh = RMF.create_rmf_file(fname)
IMP.rmf.add_hierarchy(rh, root_hier)
IMP.rmf.save_frame(rh)





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
# there are 4 restaints: Membrane inclusion for TM regions, IMS localization to localize beads within the inner radius,
# Matrix localization to keep beads farther away from outer radius and surface localization to keep the beads close to inner cylinder
# by defining max_limit and keeping the beads between the max_limit and the inner radius
# add center beads for membrane surface so that it can score based on those beads only
# a cylinder can be used (using BILD files) to visualize the restraints
#
# select_by_tuple_2(start,stop,molname,copynum,statenum);  use 'None' for them which will get all copies and states
# # # -------------------------------

############ MIC60 #################
outer_R = 50 #radius of the outer cylinder
inner_r = 30 #radius of the inner cylinder
sigma = 1 #weight
Max_limit = 10 #for surface localization, this is the distance from the center. the beads should stay within this and inner radius
thickness = outer_R-inner_r

TM_regions = (149,171,'MIC60',None,None)
mic60 = IMP.pmi.tools.select_by_tuple_2(root_hier,TM_regions,10)
mir = MembraneInclusionRestraintC(mdl,mic60, outer_R,inner_r,sigma)
output_objects.append(mir)
print("Membrane inclusion Restraint applied")

mem_surface = (627,751, 'MIC60', None, None)
mic60_19 = IMP.pmi.tools.select_by_tuple_2(root_hier, mem_surface ,10)
slr = SurfaceLocalizationRestraintC(mdl, mic60_19, inner_r, sigma, Max_limit)
output_objects.append(slr)
print("surface localization restraint applied")

mic60_C = (149,582,'MIC60',None, None)
mic60_C_ = IMP.pmi.tools.select_by_tuple_2(root_hier, mic60_C ,10)
ilr = IMSLocalizationRestraintC(mdl,mic60_C_,inner_r,sigma)
output_objects.append(ilr)
print("IMS localization restraint applied")

mic60_N = (1,148,'MIC60',None, None)
mic60_N_ = IMP.pmi.tools.select_by_tuple_2(root_hier, mic60_N ,10)
mlr = MatrixLocalizationRestraintC(mdl,mic60_N_,outer_R,sigma)
output_objects.append(mlr)
print("matrix localization restraint applied")



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

print("Excluded volume restraint applied")


# # # -----------------------------
# # # %%%%% MINIMUM PAIR RESTRAINT
# # mpr1 = IMP.pmi.restraints.basic.MinimumPairRestraint(tuple_selection1=(468,546,"MTA1",0),tuple_selection2=(1,425,"RBBP4",2),distmax=10.0,root_hier=root_hier,label="MTA1-RBBP4.2_mpr1")
# # mpr2 = IMP.pmi.restraints.basic.MinimumPairRestraint(tuple_selection1=(468,546,"MTA1",1),tuple_selection2=(1,425,"RBBP4",3),distmax=10.0,root_hier=root_hier,label="MTA1.1-RBBP4.3_mpr2")
# #
# # output_objects.append(mpr1)
# # output_objects.append(mpr2)


# # -------------------------
# # %%%%% CROSSLINKING RESTRAINT
# #
# # Restrains two particles via a distance restraint based on
# # an observed crosslink.
# #
# # First, create the crosslinking database from the input file
# # The "standard keys" correspond to a crosslink csv file of the form:
# #
# # Protein1,Residue1,Protein2,Residue2
# # A,18,G,24
# # A,18,G,146
# # A,50,G,146
# # A,50,G,171
# # A,50,G,189
# #
# # This restraint allows for ambiguity in the crosslinked residues,
# # a confidence metric for each crosslink and multiple states.
# # See the PMI documentation or the MMB book chapter for a
# # full discussion of implementing crosslinking restraints.
#
# # This first step is used to translate the crosslinking data file.
# # The KeywordsConverter maps a column label from the xl data file
# # to the value that PMI understands.
# # Here, we just use the standard keys.
# # One can define custom keywords using the syntax below.
# # For example if the Protein1 column header is "prot_1"
# # xldbkc["Protein1"]="prot_1"
#
# # The CrossLinkDataBase translates and stores the crosslink information
# # from the file "xl_data" using the KeywordsConverter.
#
# xldbkc = IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
# xldbkc.set_standard_keys()
#
# xldb = IMP.pmi.io.crosslink.CrossLinkDataBase()
# xldb.create_set_from_file(file_name=xl_data,
#                                  converter=xldbkc)
# xlr_DSSO = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
#                 root_hier=root_hier,    # Must pass the root hierarchy to the system
#                 database=xldb, # The crosslink database.
#                 length=25,              # The crosslinker plus side chain length
#                 resolution=1,           # The resolution at which to evaluate the crosslink
#                 slope=0.0001,           # This adds a linear term to the scoring function
#                 label="all_xl",                        #   to bias crosslinks towards each other
#                 weight=10)       # Scaling factor for the restraint score.
#
# output_objects.append(xlr_DSSO)
#
# xlr_BDP_PIR = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
#                 root_hier=root_hier,    # Must pass the root hierarchy to the system
#                 database=xldb, # The crosslink database.
#                 length=33,              # The crosslinker plus side chain length
#                 resolution=1,           # The resolution at which to evaluate the crosslink
#                 slope=0.0001,           # This adds a linear term to the scoring function
#                 label="all_xl",                        #   to bias crosslinks towards each other
#                 weight=10)       # Scaling factor for the restraint score.
#
# output_objects.append(xlr_BDP_PIR)
#
# print("Cross-linking restraint applied")
#
#
# #####################################################
# ###################### SAMPLING #####################
# #####################################################
# With our representation and scoring functions determined, we can now sample
# the configurations of our model with respect to the information.
print("The type of run is: " + str(runType))
print("Number of sampling frames: " + str(num_frames))


############ SHUFFLING ##################
## 3 bounding boxes for 3 regions - membrane, IMS and Matrix to shuffle the proteins in that region only

# side of the maximal square inside the cylinder

molecules = root_hier.get_children()[0].get_children()
# molecules.get_residue('MIC60',150)
# print(molecules)


IMS_span = inner_r * math.sqrt(2)
Matrix_span = outer_R * math.sqrt(2)
mem_span =  thickness *math.sqrt(2)

membrane = (mem_span/2,mem_span/2,0),(mem_span/2,mem_span/2,200)
# temp = [x for x in molecules if ('MIC60' == x.get_name())]
# print(temp)
# exit()

# shuffle_configuration(temp, bounding_box=bb_IMS, avoidcollision_rb=False)


# First shuffle all particles to randomize the starting point of the
# system. For larger systems, you may want to increase max_translation

IMP.pmi.tools.shuffle_configuration(mic60,
                                    bounding_box = membrane,
                                    avoidcollision_rb = False,
                                    max_translation=50)
                                    # excluded_rigid_bodies=fixed_set1_core)
                                    # hierarchies_included_in_collision=fixed_set1_core)

# Shuffling randomizes the bead positions. It's good to
# allow these to optimize first to relax large connectivity
# restraint scores.  100-500 steps is generally sufficient.
dof.optimize_flexible_beads(1000)



for i in output_objects:         #Add all restarints to the model
    i.add_to_model()


IMP.pmi.dof.DegreesOfFreedom.enable_all_movers(dof)

# Now, add all of the other restraints to the scoring function to start sampling
evr.add_to_model()
# xlr.add_to_model()

print("Replica Exchange Maximum Temperature : " + str(rex_max_temp))

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange0(mdl,
        root_hier=root_hier,                    # pass the root hierarchy
        # crosslink_restraints=[xlr],
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
