from af_pipeline.RigidBodies import RigidBodies
import sys
from af_pipeline.Interaction import Interaction

structure_path = sys.argv[1]
data_path = sys.argv[2]
af_rigid = RigidBodies(
    structure_path=structure_path,
    data_path=data_path,
    average_atom_pae=True
   
)

# parameters to vary
af_rigid.plddt_cutoff = 70
af_rigid.plddt_cutoff_idr = 50 # you can set different cutoff for IDR
af_rigid.pae_cutoff = 12 # Edges will be created between all the residues < PAE cutoff
af_rigid.pae_power = 1
af_rigid.resolution = 0.5 # default value in ChimeraX
# lower value of resolution results in larger domains
af_rigid.library = "igraph" # "networkx" is slower

domains = af_rigid.predict_domains(
    num_res=5, # minimum number of residues in a domain
    num_proteins=2, # minimum number of proteins in a domain
    plddt_filter=True, # filter domains based on pLDDT score
)

af_rigid.save_rigid_bodies(
    domains=domains,
    output_dir="./outputs/",
    output_format="txt",
    save_structure=True, # if set to True, you will get each rigid body as a separate PDB
    no_plddt_filter_for_structure=False, # if set to True, pLDDT filter will be ignored for saving the PDB
)

# average pLDDT
af_rigid.get_average_pLDDT(domains=domains, chain_type="any")
# valid values for chain_type are: "any", "idr", "r"

all_interface_residues = af_rigid.get_interface_residues(
    domains=domains,
    contact_threshold=8,
    as_matrix=False
)

af_rigid.get_ipLDDT(
  all_interface_residues=all_interface_residues,
  interface_type="any-any"
)

all_interface_residues = af_rigid.get_interface_residues(
    domains=domains,
    contact_threshold=8,
    as_matrix=True
)

af_rigid.get_ipae(
    all_interface_residues=all_interface_residues
)

regions_of_interest_ = []

af_interaction = Interaction(
    struct_file_path=structure_path,
    data_file_path=data_path,
    output_dir='./outputs/',
    average_atom_pae=True
)

af_interaction.plddt_cutoff = 70
af_interaction.pae_cutoff = 5
af_interaction.interaction_map_type = "contact"
af_interaction.contact_threshold = 8

if not regions_of_interest_:
    regions_of_interest_ = af_interaction.create_regions_of_interest()

print(f"Regions of interest for {structure_path}:")
for i, region_of_interest in enumerate(regions_of_interest_):
    print(f"Region {i + 1}: {region_of_interest}")
print("-------------------------------")

for region_of_interest in regions_of_interest_:
    af_interaction.save_ppair_interaction(
        region_of_interest=region_of_interest,
        save_plot=True,
        plot_type="both",
        concat_residues=True,
        contact_probability=True,
    )
print("-------------------------------")

