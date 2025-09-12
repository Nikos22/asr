This repository contains scripts and links to data used in Vakirlis et al. 2024, https://doi.org/10.1093/gbe/evae151

If you use any of these scripts, please cite: 
Vakirlis, N., Acar, O., Cherupally, V. & Carvunis, A.-R. Ancestral Sequence Reconstruction as a Tool to Detect and Study De Novo Gene Emergence. Genome Biology and Evolution 16, evae151 (2024).


The script run_reconstructions_clean.sh runs the necessary tools to obtain the reconstructions.
You run it like this: ./run_reconstructions_clean.sh YBR296C-A.nt.fa species_tree.nwk 8 
in which argument one contains the input sequences in FASTA format, argument two is the species tree in newick format, and argument three is the number of species in your analysis (and also in your species tree).
This is needed because if a given input sequences is missing some species, the script adjusts. Please see comments etc. within the script. You also need to specify your outgroup species (see line 8).

To run the script you need the following:
- a parser called Fast2Phylip.pl (found here)
- MAFFT
- raxml-ng
- Gotree and Goalign
- PhyloFit and PREQUEL from PHAST package
- PRANK
- FastML (you have to set the path for that one in the script, see lines 71 and 72)

The output of all the tools is then found in a directory with your input sequence name (e.g. here would be YBR296C-A).

The second script is called anc_RFC.py and what it does is that it goes ancestor by ancestor for the various reconstructions (output of previous script), and finds the ancestral ORF with the best RFC score relative to the ORF in the focal species. 
You run it like this: python anc_RFC.py [DIR_NAME] [FOCAL_SPECIES] [NO_OF_RANDOMIZATIONS]
The first argument is the directory that the previous script creates, where all the files are stored.
The second one is the name of the focal species, as it appears within the initial FASTA file (e.g. here would be Scer).
The third is the number of randomizations that the script will perform to calculate the empirical P-value.
The script outputs a table with the RFC and P-value for each ancestor and tool, and it also stores all the best RFC ORFs in a file with the suffix "_bestRFC_ORFs.fasta" .

The main output looks like the following:

###
YBR296C-A FastML_joint species N1 0 Sarb_Scer_Seub_Sjur_Skud_Smik_Spar_Suva 0.35 0.0
YBR296C-A FastML_joint species N2 1 Sarb_Scer_Sjur_Skud_Smik_Spar 0.35 0.0
YBR296C-A FastML_joint species N3 2 Scer_Sjur_Skud_Smik_Spar 0.717 0.0
YBR296C-A FastML_joint species N4 3 Scer_Sjur_Smik_Spar 0.717 0.0
YBR296C-A FastML_joint species N5 4 Scer_Sjur_Spar 0.717 0.0
YBR296C-A FastML_joint species N6 5 Scer_Spar 0.717 0.0
YBR296C-A FastML_joint species Scer_spec NA Scer -1 -1
###
(disregard the last line)

The column names are the following:
"gene_name", "tool", "topology", "ancestor_code_of_tool", "ancestor_order", "ancestor_descendants", "best_RFC", "P-value"
