# Smc-Hawks

Scripts used to generate network data from raw hhblits/hhsearch output files. Some additional bits and pieces, mostly for dealing with fasta data and parsing other HHsuite output files. The bulk of the work is carried out by hhr_network.jl and control_networks.jl. The former generates the networks from the results of hhsuite searches, and carries out all of the quality control and normalisation steps. The latter does the computationally intensive permutation tests for the network clusters.

supp_data file contains the following:
* network_seeds.txt - Uniprot IDs of the sequences used to initialise network searches (see paper supplemental methods for details).
* species_network.txt - Processed network files generated from raw hhblits output files.
* mcl_clusters.txt - List of the cluster IDs produced from running MCL algorithm on the networks.
* go_cluster_enrichments.csv - GO biological process cluster enrichment data from a selection of human network clusters.
* eukaryote_species_tested.txt - A selection of eukaryote species in which high confidence hawk proteins were found. less confident ones can be found by searching eggnog.
* software_used.docx - table of all the software and parameters used in the work.

If you would like access to any of the raw data then drop me a message and I will be happy to help. Likewise if you have any questions about the scripts.
