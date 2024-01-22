
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
run_output_dir = 'run_' + runID
data_direc = sys.argv[3]

if runType == "test":
    num_frames = 1500
elif runType == "prod":
    num_frames = 30000

# Topology File
topology_file = f'{data_direc}/topology_mic19_1full_1N_1C_independent.txt'

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
                                  max_srb_trans= 0.4,
                                  max_srb_rot=0.02)

molecules = t.get_components()



rex_max_temp = 1.25

# mem_surface = [(634, 637,'MIC60', None, None),(643, 644,'MIC60', None, None),(693,693,'MIC60', None, None)]
# fixed_particles = IMP.pmi.tools.select_by_tuple_2(root_hier,(627,758,'MIC60', None, None),1)
fixed_particles = []


fixed_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC60',copy_indexes = [0,1],residue_indexes = range(627,753)).get_selected_particles()
fixed_particles += IMP.atom.Selection(hierarchy = root_hier,molecule='MIC19',copy_indexes= [0,2],residue_indexes = range(186,228)).get_selected_particles()
# doesnot work if select residue ranges, can work with full protein


print(fixed_particles)
IMP.pmi.dof.DegreesOfFreedom.enable_all_movers(dof)

fixed_beads,fixed_rbs=dof.disable_movers(fixed_particles, [IMP.core.RigidBodyMover,
                                         IMP.pmi.TransformMover])

print(fixed_rbs)
# exit()

IMP.pmi.tools.shuffle_configuration(root_hier,
                                    max_translation=50,
                                    excluded_rigid_bodies=fixed_rbs)

dof.optimize_flexible_beads(1000)


print("Replica Exchange Maximum Temperature : " + str(rex_max_temp))

# Run replica exchange Monte Carlo sampling
rex=IMP.pmi.macros.ReplicaExchange(mdl,
        root_hier=root_hier,                    # pass the root hierarchy
        # crosslink_restraints=[xlr],
        # This allows viewing the crosslinks in Chimera. Also, there is not inter-protein ADH crosslink available. Hence it is not mentioned in this list
        monte_carlo_temperature = 1.0,
        replica_exchange_minimum_temperature = 1.0,
        replica_exchange_maximum_temperature = rex_max_temp,
	    monte_carlo_sample_objects=dof.get_movers(),  # pass all objects to be moved ( almost always dof.get_movers() )
        global_output_directory=run_output_dir,      # The output directory for this sampling run.
        # output_objects=output_objects,          # Items in output_objects write information to the stat file.
        monte_carlo_steps=10,                   # Number of MC steps between writing frames
        number_of_best_scoring_models=0,        # set >0 to store best PDB files (but this is slow)
        number_of_frames=num_frames)            # Total number of frames to run / write to the RMF file.
        #test_mode=test_mode)                    # (Ignore this) Run in test mode (don't write anything)

# Ok, now we finally do the sampling!
rex.execute_macro()
#
# Outputs are then analyzed in a separate analysis script.
