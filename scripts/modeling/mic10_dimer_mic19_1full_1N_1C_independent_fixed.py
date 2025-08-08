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
import pandas as pd

# ### --------------------------------
### WRAPPERS ###
### --------------------------------


# Wrapper for the Biochemical Binding Restraint
class MinimumPairDistanceBindingRestraint(IMP.pmi.restraints.RestraintBase):
    def __init__(self, model, plist1, plist2, x0=0, kappa=1, label=None, weight=1):
        name = "MinimumPairDistanceBindingRestraint%1%"
        super(MinimumPairDistanceBindingRestraint, self).__init__(
            model, name=name, label=label, weight=weight
        )
        l1 = IMP.container.ListSingletonContainer(mdl)
        l1.add(plist1)
        l2 = IMP.container.ListSingletonContainer(mdl)
        l2.add(plist2)
        bipartite_container = IMP.container.AllBipartitePairContainer(l1, l2)
        score = IMP.core.HarmonicSphereDistancePairScore(x0, kappa)
        res_main = IMP.container.MinimumPairRestraint(score, bipartite_container, 1)
        self.rs.add_restraint(res_main)
        print(
            "RESTRAINT: Added MPDBR on particles",
            len(plist1),
            ":",
            len(plist2),
            "at x0:kappa",
            str(x0),
            ":",
            str(kappa),
        )


# Note: sigma is mem_wt
# Wrapper for Transmembrane Restraint
class TransMembraneRestraintC(IMP.pmi.restraints.RestraintBase):
    def __init__(self, model, plist, R, r, sigma, weight=1):
        particles = plist
        name = "TransMembraneRestraint%1%"
        super(TransMembraneRestraintC, self).__init__(model, name=name, weight=weight)
        res_main = IMP.micos.TransMembraneRestraint(particles, R, r, sigma)
        self.rs.add_restraint(res_main)

# Wrapper for IMS localization restraint
class IMSLocalizationRestraintC(IMP.pmi.restraints.RestraintBase):
    def __init__(self, model, plist, r, sigma, weight=1):
        particles = plist
        name = "IMSLocalizationRestraint%1%"
        super(IMSLocalizationRestraintC, self).__init__(model, name=name, weight=weight)
        res_main = IMP.micos.IMSLocalizationRestraint(particles, r, sigma)
        self.rs.add_restraint(res_main)


# Wrapper for matrix localization restraint
class MatrixLocalizationRestraintC(IMP.pmi.restraints.RestraintBase):
    def __init__(self, model, plist, R, sigma, weight=1):
        particles = plist
        name = "MatrixLocalizationRestraint%1%"
        super(MatrixLocalizationRestraintC, self).__init__(
            model, name=name, weight=weight
        )
        res_main = IMP.micos.MatrixLocalizationRestraint(particles, R, sigma)
        self.rs.add_restraint(res_main)


# Wrapper for zaxial restraint on protein level
class ZAxialRestraintC(IMP.pmi.restraints.RestraintBase):
    def __init__(self, model, plist, ub, lb, sigma, method, label=None, weight=1):
        particles = plist
        name = "ZAxialRestraint%1%"
        super(ZAxialRestraintC, self).__init__(
            model, name=name, weight=weight, label=label
        )
        res_main = IMP.micos.ZAxialRestraint(particles, ub, lb, sigma, method)
        self.rs.add_restraint(res_main)


runType = sys.argv[1]  # Specify test or prod
runID = sys.argv[2]  # Specify the number of runs
run_output_dir = "run_" + runID
data_direc = sys.argv[3]

if runType == "test":
    num_frames = 2000
elif runType == "prod":
    num_frames = 50000

# Xlinkdata files
xl_BDP_PIR_human = f"{data_direc}/crosslinks/human/sampling_BDP_PIR_human.csv"
xl_BDP_PIR_mouse = f"{data_direc}/crosslinks/human/sampling_BDP_PIR_mouse.csv"
xl_DSS_BS3_human = f"{data_direc}/crosslinks/human/sampling_DSS_BS3_human.csv"
xl_DHSO_DSSO_human = f"{data_direc}/crosslinks/human/sampling_DHSO_DSSO_human.csv"  # yu. bartolec, xlinkdb dsso xlinks
xl_DSSO_mouse = f"{data_direc}/crosslinks/human/sampling_DSSO_mouse.csv"

# AF3 predictions and biochemical data file
af3_predictions_biochemical = f"{data_direc}/biochemical/af3_predictions_biochemical.csv"

# Topology File
topology_file = f"{data_direc}/topology_mic10_dimer_mic19_1full_1N_1C_independent.txt"

# Weights
mem_wt = 0.04
zar_wt = 0.2
coIP_wt = 1
AF_wt = 2


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
root_hier, dof = bs.execute_macro(max_rb_trans= 0.5,
                                  max_rb_rot= 0.1,
                                  max_bead_trans= 2,
                                  max_srb_trans= 0.4,
                                  max_srb_rot=0.02)

rex_max_temp = 1.67

########### Fixing Mic60-Mic19 tetramer close to the membrane ######################
fixed_particles = []


fixed_particles += IMP.atom.Selection(
    hierarchy=root_hier,
    molecule="MIC60",
    copy_indexes=[0, 1],
    residue_indexes=range(627, 752),
).get_selected_particles()
fixed_particles += IMP.atom.Selection(
    hierarchy=root_hier, molecule="MIC19", residue_indexes=range(186, 227)
).get_selected_particles()
# doesnot work if select residue ranges, can work with full protein
# print(fixed_particles)


# Uncomment the following lines to get test.rmf file to visualise the system representation

# Uncomment this line for verbose output of the representation
# IMP.atom.show_with_representations(root_hier)
# # output to RMF
# fname = 'test.rmf'
# rh = RMF.create_rmf_file(fname)
# IMP.rmf.add_hierarchy(rh, root_hier)
# IMP.rmf.save_frame(rh)

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
# Matrix localization to keep beads farther away from outer radius
# a cylinder can be used (using BILD files) to visualize the restraints
# # # -------------------------------

## These selections are for the membrane restraints.
protien_selections = {'mic10_n': IMP.atom.Selection(hierarchy=root_hier, molecule="MIC10", copy_indexes=[0, 1], resolution=10,
                                                    residue_indexes=range(1, 13),).get_selected_particles(),

    'mic10_tm1': IMP.atom.Selection(hierarchy=root_hier, molecule="MIC10", copy_indexes=[0, 1], resolution=10,
                                    residue_indexes=range(13, 37),).get_selected_particles(),

    'mic10_matrix': IMP.atom.Selection(hierarchy=root_hier, molecule="MIC10", copy_indexes=[0, 1], resolution=10,
                                      residue_indexes=range(37, 40),).get_selected_particles(),

    'mic10_tm2': IMP.atom.Selection(hierarchy=root_hier, molecule="MIC10", copy_indexes=[0, 1], resolution=10,
                                    residue_indexes=range(40, 61),).get_selected_particles(),

    'mic10_c': IMP.atom.Selection(hierarchy=root_hier, molecule="MIC10", copy_indexes=[0, 1], resolution=10,
                                  residue_indexes=range(61, 79),).get_selected_particles(), 
    
    'mic13_n': IMP.atom.Selection(hierarchy=root_hier, molecule="MIC13", resolution=10, 
                                  residue_indexes=range(1, 8),).get_selected_particles(),
    
    'mic13_tm': IMP.atom.Selection(hierarchy=root_hier, molecule="MIC13", resolution=10, 
                                   residue_indexes=range(8, 24),).get_selected_particles(),

    'mic13_central_c': IMP.atom.Selection(hierarchy=root_hier, molecule="MIC13", resolution=10,
                                          residue_indexes=range(24, 119),).get_selected_particles(),

    'mic60_cc': IMP.atom.Selection(hierarchy=root_hier, molecule="MIC60", copy_indexes=[0, 1, 2, 3], resolution=10, 
                                   residue_indexes=range(410, 583),).get_selected_particles(),

    'mic60_link': IMP.atom.Selection(hierarchy=root_hier, molecule="MIC60", copy_indexes=[0, 1], resolution=10, 
                                    residue_indexes=range(583, 627),).get_selected_particles(),

    'mic60_lbs_m': IMP.atom.Selection(hierarchy=root_hier, molecule="MIC60", copy_indexes=[0, 1], resolution=10, 
                                      residue_indexes=range(627, 759),).get_selected_particles(),

    'mic19_n_cc': IMP.atom.Selection(hierarchy=root_hier, molecule="MIC19", copy_indexes=[0, 1], resolution=10, 
                                     residue_indexes=range(1, 186),).get_selected_particles(),

    'mic19_chch': IMP.atom.Selection(hierarchy=root_hier, molecule="MIC19", copy_indexes=[0, 2], resolution=10,
                                     residue_indexes=range(186, 228),).get_selected_particles(),

    'mic19_n_sam50': IMP.atom.Selection(hierarchy=root_hier, molecule="MIC19", copy_indexes=[0, 1], resolution=10,
                                        residue_indexes=range(1, 15),).get_selected_particles()}

R = 130  # radius of the outer cylinder
r = 90  # radius of the inner cylinder
on_cj = 0
max_in_cj = 50
max_above_cj = -30

# TransMembraneRestraint
for particle_set in (protien_selections['mic10_tm1'], protien_selections['mic10_tm2'], protien_selections['mic13_tm']):
    tmr = TransMembraneRestraintC(mdl, particle_set, R, r, mem_wt)
    output_objects.append(tmr)
    tmr.add_to_model()

# IMSLocalizationRestraint
for particle_set in (
    protien_selections['mic10_n'],
    protien_selections['mic10_c'],
    protien_selections['mic13_central_c'],
    protien_selections['mic60_cc'],
    protien_selections['mic60_link'],
    protien_selections['mic60_lbs_m'],
    protien_selections['mic19_n_cc'],
    protien_selections['mic19_chch'],
):
    ilr = IMSLocalizationRestraintC(mdl, particle_set, r, mem_wt)
    output_objects.append(ilr)
    ilr.add_to_model()

# MatrixLocalizationRestraint
for particle_set in (protien_selections['mic10_matrix'], protien_selections['mic13_n'],):
    mlr = MatrixLocalizationRestraintC(mdl, particle_set, R, mem_wt)
    output_objects.append(mlr)
    mlr.add_to_model()

# ZAxialRestraint for interaction with SAM/TOB complex
## for Mic60 ateast one bead above cj
for particle_set in (protien_selections['mic19_n_sam50'], 
                     protien_selections['mic60_cc'], 
                     protien_selections['mic60_link'],
                     protien_selections['mic60_lbs_m']):
    zar_sam50 = ZAxialRestraintC(
        mdl, particle_set, max_above_cj, max_in_cj, zar_wt, "domain", label="zar_sam50"
    )
    output_objects.append(zar_sam50)
    zar_sam50.add_to_model()


# ZAxialRestraint in the cylinder
for particle_set in (
    protien_selections['mic10_n'],
    protien_selections['mic10_tm1'],
    protien_selections['mic10_matrix'],
    protien_selections['mic10_tm2'],
    protien_selections['mic10_c'],
    protien_selections['mic13_n'],
    protien_selections['mic13_tm'],
    protien_selections['mic13_central_c'],
    
):
    zar_below = ZAxialRestraintC(
        mdl, particle_set, on_cj, max_in_cj, zar_wt, "average", label="zar_below"
    )
    output_objects.append(zar_below)
    zar_below.add_to_model()


# To specify the particles in the specific bounding box
tm_particles = []
matrix_particles = []
ims_mem_surface_particles = []

ims_mem_surface_particles.extend(
    [
        protien_selections['mic10_n'],
        protien_selections['mic10_c'],
        protien_selections['mic13_central_c'],
        protien_selections['mic60_cc'], 
        protien_selections['mic60_link'],
        protien_selections['mic60_lbs_m'],
        protien_selections['mic19_n_cc'],
        protien_selections['mic19_chch'],
    ]
)
tm_particles.extend([protien_selections['mic10_tm1'], protien_selections['mic10_tm2'], protien_selections['mic13_tm'],])
matrix_particles.extend([protien_selections['mic10_matrix'], protien_selections['mic13_n'],])

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
    included_objects=[root_hier], resolution=1000
)
output_objects.append(evr)
print("Excluded volume restraint applied")


# #####################################################
# ###################### SAMPLING #####################
# #####################################################
# With our representation and scoring functions determined, we can now sample
# the configurations of our model with respect to the information.

############ SHUFFLING ##################
## 3 bounding boxes for 3 regions - membrane, IMS and Matrix to shuffle the proteins in that region only
ims_bb = ((0, -r, 0), (r, r, 100))  # this is for half cylinder
tm_bb = ((r, -r, 0), (R, R, 100))

# First shuffle all particles to randomize the starting point of the
# system. For larger systems, you may want to increase max_translation
IMP.pmi.dof.DegreesOfFreedom.enable_all_movers(dof)
fixed_beads, fixed_rbs = dof.disable_movers(
    fixed_particles, [IMP.core.RigidBodyMover, IMP.pmi.TransformMover]
)
IMP.pmi.tools.shuffle_configuration(
    ims_mem_surface_particles,
    max_translation=100,
    bounding_box=ims_bb,
    excluded_rigid_bodies=fixed_rbs,
    hierarchies_included_in_collision=fixed_particles,
)

IMP.pmi.tools.shuffle_configuration(
    tm_particles, max_translation=100, bounding_box=tm_bb, avoidcollision_rb=False
)

IMP.pmi.tools.shuffle_configuration(matrix_particles, max_translation=100)
# Shuffling randomizes the bead positions. It's good to
# allow these to optimize first to relax large connectivity
# restraint scores.  100-500 steps is generally sufficient.
dof.optimize_flexible_beads(1000)


print("Replica Exchange Maximum Temperature : " + str(rex_max_temp))


# -----------------------------
# %%%%% MINIMUM PAIR RESTRAINT

df = pd.read_csv(af3_predictions_biochemical)

for idx, row in df.iterrows():
    p1 = row['Protein1']
    p2 = row['Protein2']
    r1 = row['Region1']
    r2 = row['Region2']
    
    wt = coIP_wt if idx == 0 else AF_wt

    if '-' in r1:
        start1, end1 = map(int, r1.split('-'))
    else:
        start1 = end1 = int(r1)
    
    if '-' in r2:
        start2, end2 = map(int, r2.split('-'))
    else:
        start2 = end2 = int(r2)

    mpr = MinimumPairDistanceBindingRestraint(
    mdl,
    IMP.pmi.tools.select_by_tuple_2(root_hier, (start1, end1, p1, None, None), 1),
    IMP.pmi.tools.select_by_tuple_2(root_hier, (start2, end2, p2, None, None), 1),
    5,
    1,
    f'{p1}_{p2}_mpr',
    wt, )
    output_objects.append(mpr)
    mpr.add_to_model()

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
xldb_BDP_PIR_human.create_set_from_file(file_name=xl_BDP_PIR_human, converter=xldbkc)
xlr_BDP_PIR_human = (
    IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
        root_hier=root_hier,  # Must pass the root hierarchy to the system
        database=xldb_BDP_PIR_human,  # The crosslink database.
        length=42,  # The crosslinker plus side chain length
        resolution=1,  # The resolution at which to evaluate the crosslink
        slope=0.0001,  # This adds a linear term to the scoring function
        label="BDP_PIR_human",  #   to bias crosslinks towards each other
        weight=5,  # Scaling factor for the restraint score.
        linker=ihm.ChemDescriptor("bruce"),
    )
)

output_objects.append(xlr_BDP_PIR_human)


xldb_BDP_PIR_mouse = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_BDP_PIR_mouse.create_set_from_file(file_name=xl_BDP_PIR_mouse, converter=xldbkc)
xlr_BDP_PIR_mouse = (
    IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
        root_hier=root_hier,  # Must pass the root hierarchy to the system
        database=xldb_BDP_PIR_mouse,  # The crosslink database.
        length=42,  # The crosslinker plus side chain length
        resolution=1,  # The resolution at which to evaluate the crosslink
        slope=0.0001,  # This adds a linear term to the scoring function
        label="BDP_PIR_mouse",  #   to bias crosslinks towards each other
        weight=4,  # Scaling factor for the restraint score.
        linker=ihm.ChemDescriptor("bruce"),
    )
)

output_objects.append(xlr_BDP_PIR_mouse)


xldb_DSS_BS3 = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_DSS_BS3.create_set_from_file(file_name=xl_DSS_BS3_human, converter=xldbkc)
xlr_DSS_BS3 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=root_hier,  # Must pass the root hierarchy to the system
    database=xldb_DSS_BS3,  # The crosslink database.
    length=25,  # The crosslinker plus side chain length
    resolution=1,  # The resolution at which to evaluate the crosslink
    slope=0.0001,  # This adds a linear term to the scoring function
    label="DSS_BS3",  #   to bias crosslinks towards each other
    weight=5,  # Scaling factor for the restraint score.
    linker=ihm.ChemDescriptor("bruce"),
)

output_objects.append(xlr_DSS_BS3)


xldb_DHSO_DSSO = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_DHSO_DSSO.create_set_from_file(file_name=xl_DHSO_DSSO_human, converter=xldbkc)
xlr_DHSO_DSSO = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=root_hier,  # Must pass the root hierarchy to the system
    database=xldb_DHSO_DSSO,  # The crosslink database.
    length=21,  # The crosslinker plus side chain length
    resolution=1,  # The resolution at which to evaluate the crosslink
    slope=0.0001,  # This adds a linear term to the scoring function
    label="DHSO_DSSO",  #   to bias crosslinks towards each other
    weight=5,  # Scaling factor for the restraint score.
    linker=ihm.ChemDescriptor("bruce"),
)

output_objects.append(xlr_DHSO_DSSO)


xldb_DSSO_mouse = IMP.pmi.io.crosslink.CrossLinkDataBase()
xldb_DSSO_mouse.create_set_from_file(file_name=xl_DSSO_mouse, converter=xldbkc)
xlr_DSSO_mouse = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(
    root_hier=root_hier,  # Must pass the root hierarchy to the system
    database=xldb_DSSO_mouse,  # The crosslink database.
    length=21,  # The crosslinker plus side chain length
    resolution=1,  # The resolution at which to evaluate the crosslink
    slope=0.0001,  # This adds a linear term to the scoring function
    label="DSSO_mouse",  #   to bias crosslinks towards each other
    weight=4,  # Scaling factor for the restraint score.
    linker=ihm.ChemDescriptor("bruce"),
)
output_objects.append(xlr_DSSO_mouse)

xlr_BDP_PIR_human.add_to_model()
xlr_BDP_PIR_mouse.add_to_model()
xlr_DSS_BS3.add_to_model()
xlr_DHSO_DSSO.add_to_model()
xlr_DSSO_mouse.add_to_model()

evr.add_to_model()
# Run replica exchange Monte Carlo sampling
rex = IMP.pmi.macros.ReplicaExchange(
    mdl,
    root_hier=root_hier,  # pass the root hierarchy
    monte_carlo_temperature=1.0,
    replica_exchange_minimum_temperature=1.0,
    replica_exchange_maximum_temperature=rex_max_temp,
    monte_carlo_sample_objects=dof.get_movers(),  # pass all objects to be moved ( almost always dof.get_movers() )
    global_output_directory=run_output_dir,  # The output directory for this sampling run.
    output_objects=output_objects,  # Items in output_objects write information to the stat file.
    monte_carlo_steps=10,  # Number of MC steps between writing frames
    number_of_best_scoring_models=0,  # set >0 to store best PDB files (but this is slow)
    number_of_frames=num_frames,
)

rex.execute_macro()