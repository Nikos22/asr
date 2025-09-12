#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  anc_RFC.py
#  
#  Copyright 2025 Nikolaos Vakirlis
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
# 
# to run the script:
# python anc_RFC.py [DIR_NAME] [FOCAL_SPECIES] [NO_OF_RANDOMIZATIONS] 
#

##IMPORT 

import string, sys,os, re, os, random, copy
import numpy as np
from Bio import SeqIO, Phylo
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices
from os import listdir
from os.path import join , isfile, basename


global subs_mat
subs_mat = substitution_matrices.load('TRANS')


def get_best_orf(alignment, scer_seq, anc_seq, gene_name) :
		#this function maps ORFs of an ancestral sequence on to the main ORF (scer_seq) based on a pairwise alignment of the two
		#it returns the ORF with the best reading frame conservation
		
		frame_cod_scer = []
		frame_cod_anc = []
		frame_scer=0
		frame_anc=0
		cfs_by_length = []
		ATG_pos_scer = None
		gap_count_anc = 0
		ATG_pos_anc = []
		stop_pos_anc = []
		
		best_anc_ORF = None
		anc_ORFs = []
		anc_ORFs_tr_co = []
		
		anc_ORF_start_pos = 0
		anc_ORF_stop_pos = 0
		for p in range(0, len(alignment[0][0])) :
			
			
			if alignment[0][0][p] == '-' :
				frame_scer=(frame_scer+1)%3
				frame_cod_scer.append("-")
				gap_count_anc+=1
			else :
				frame_cod_scer.append(frame_scer)
			
			if alignment[0][0][p:p+3] == 'ATG' and ATG_pos_scer==None :
				ATG_pos_scer=p
				
			if alignment[0][1][p] == '-' :
				frame_anc=(frame_anc+1)%3
				frame_cod_anc.append("-")
			else :
				frame_cod_anc.append(frame_anc)
		
		SeqIO.write( anc_seq, "./gene.temp.fasta", "fasta")
		# to use definition of ORFs START-STOP replace "-find 2" with "-find 3"
		os.system("getorf -find 2 -minsize 6 -sequence ./gene.temp.fasta -outseq ./gene.temp.orfs 2> /dev/null")
		os.system("grep '>' ./gene.temp.orfs | grep -v REVERSE |cut -f 2 -d '[' | tr -d ']-' > ./gene.temp.orfs.clean")
		
		with open("./gene.temp.orfs.clean") as f :
			for l in f :
				seq_start = int(l.split()[0])-1
				seq_end = int(l.split()[1])+3
				counter_start = 0
				counter_end = 0
				found_start = False
				found_end = False
				for p in range(len(alignment[0][1])) :
					
						
					if alignment[0][1][p] != '-' :
						if counter_start<seq_start :
							counter_start+=1
						elif not found_start :
							counter_start=p
							found_start=True
						if counter_end<seq_end :
							counter_end+=1
						elif not found_end :
							counter_end=p
							found_end=True	
							
					if found_start and found_end :
						break
						
				anc_ORFs.append((counter_start, counter_end))
				anc_ORFs_tr_co.append((seq_start, seq_end))
		
		max_cfs = 0
		max_cfs_orf_len = 0
		max_cfs_coords = ()
		max_orf_ord = -1
		for o in range(len(anc_ORFs)) :
			orf = anc_ORFs[o]
			ori_orf = anc_ORFs_tr_co[o]
			cfs = 0
			for p in range(orf[0], orf[1]) :
				if frame_cod_anc[p] == frame_cod_scer[p] and frame_cod_anc[p] != '-' :
					cfs+=1
			if cfs>max_cfs :
				max_cfs = cfs
				max_cfs_coords = (ori_orf[0], ori_orf[1])
				max_cfs_orf_len = ori_orf[1]-ori_orf[0]
				max_orf_ord = o
			elif cfs==max_cfs :
				if (ori_orf[1]-ori_orf[0]) > max_cfs_orf_len :
					max_cfs_orf_len = ori_orf[1]-ori_orf[0]
					max_cfs = cfs
					max_cfs_coords = (ori_orf[0], ori_orf[1])
					max_orf_ord = o
	
		#print(gene_name, anc_seq.id, max_cfs_orf_len, max_cfs, float(max_cfs/len(scer_seq)))
		
		orf_counter=0
		best_orf_seq = ""
		for Seq in SeqIO.parse(open("./gene.temp.orfs"), 'fasta') :
			if orf_counter == max_orf_ord :
				best_orf_seq=Seq
				break
			else :
				orf_counter+=1
		
		return float(round(max_cfs/len(scer_seq),3)), best_orf_seq


def randomize_and_getcfs(scer_seq, anc_seq, gene_name, iterations) :
	# this function performes the randomizations runs get_best_orf() 
	rand_anc_seq = copy.deepcopy(anc_seq)
	max_orf_list = []
	for i in range(0, iterations) :
		rand_seq = list(anc_seq.seq._data)
		random.shuffle(rand_seq)
		rand_anc_seq.seq._data=''.join(rand_seq)
		alignment = pairwise2.align.globalds(scer_seq, rand_anc_seq.seq._data, subs_mat, -3, -.1, one_alignment_only=True)
		res = get_best_orf(alignment, scer_seq, rand_anc_seq, gene_name)
		max_orf_list.append(res[0])

	return(max_orf_list)



def main(args):
	
	# this function parses and analyzes the output of run_reconstructions.sh for each ORF/gene
	# this output should be stored in a folder that has the name of the ORF/gene
	
	# the name of the ORF/gene that is also the base directory
	gene_name = args[1]
	
	# input alignment file
	ali = join(args[1], f"{args[1]}_ali.fa")
	
	# the sequence of the ORF in teh focal species
	scer_seq = None
	
	
	comp_scer_acc = None
	
	# Focal species name
	focal_name = args[2]

	# Check if the FASTA file exists before parsing
	if isfile(ali):
	    for Seq in SeqIO.parse(open(ali), 'fasta'):
	        # Convert sequence to uppercase regardless of its original case
	        Seq.seq = Seq.seq.upper()
	        
	        if Seq.id == focal_name: # this should be changed to the name of the focal species as it appears in the alignment
	            scer_seq = ''.join([x for x in Seq.seq if x != '-'])
	            Seq.seq._data = ''.join([x for x in Seq.seq if x != '-'])
	            comp_scer_acc = Seq
	            break
	else:
	    print(f"Warning: The file {ali} does not exist. Skipping sequence processing.")

	# Number of randomized iterations performed
	no_reps = int(args[3])

	# Define paths for various files
	free_tree = join(args[1], f"{args[1]}_freeTop.raxml.bestTree_midRoot")
	sp_tree = join(args[1], f"{args[1]}_fixedTop.raxml.bestTree_CORPAR")

	free_tree_for_preq = join(args[1], "temp.tree.midroot")
	sp_tree_for_preq = join(args[1], "temp.tree.corout")

	prequel_free = join(args[1], "prequel_freeTop/")
	prequel_species = join(args[1], "prequel_spTop/")

	# Check for existence of prequel files
	if os.path.isdir(prequel_free):
	    prequel_free_names = [f.split('.')[1] for f in os.listdir(prequel_free) if f.endswith(".fa")]
	else:
	    print(f"Warning: The directory {prequel_free} does not exist.")

	if os.path.isdir(prequel_species):
	    prequel_sp_names = [f.split('.')[1] for f in os.listdir(prequel_species) if f.endswith(".fa")]
	else:
	    print(f"Warning: The directory {prequel_species} does not exist.")

	prank_free = join(args[1], "prank_freeTop/")
	prank_species = join(args[1], "prank_spTop/")

	# Check and read prank ancestral trees
	if isfile(join(prank_free, "ORF_alignment.anc.dnd")):
	    prank_anc_tree_free = Phylo.read(open(join(prank_free, "ORF_alignment.anc.dnd")), "newick")
	else:
	    print(f"Warning: The file {join(prank_free, 'ORF_alignment.anc.dnd')} does not exist.")

	if isfile(join(prank_species, "ORF_alignment.anc.dnd")):
	    prank_anc_tree_sp = Phylo.read(open(join(prank_species, "ORF_alignment.anc.dnd")), "newick")
	else:
	    print(f"Warning: The file {join(prank_species, 'ORF_alignment.anc.dnd')} does not exist.")

	# Define paths for prank ancestral sequences
	prank_anc_seqs_free = join(prank_free, "ORF_alignment.anc.fas")
	prank_anc_seqs_sp = join(prank_species, "ORF_alignment.anc.fas")

	fastml_free = join(args[1], "fastml_freeTop/")
	fastml_species = join(args[1], "fastml_spTop/")

	# Handling fastml files with checks for existence
	if isfile(join(fastml_free, "seq.marginal_IndelAndChars.txt")):
	    fastml_anc_seqs_marg_free = join(fastml_free, "seq.marginal_IndelAndChars.txt") 
	    fastml_pp_free = join("~/asr/FINAL/post_probs/fastml_free/", gene_name + "_fastml_free_char_indel.txt")
	else:
	    fastml_anc_seqs_marg_free = join(fastml_free, "seq.marginal.txt")
	    fastml_pp_free = join("~/asr/FINAL/post_probs/fastml_free/", gene_name + "_fastml_free_char_only.txt")

	fastml_anc_seqs_joint_free = join(fastml_free, "seq.joint.txt")

	if isfile(join(fastml_species, "seq.marginal_IndelAndChars.txt")):
	    fastml_anc_seqs_marg_sp = join(fastml_species, "seq.marginal_IndelAndChars.txt") 
	    fastml_pp_sp = join("~/asr/FINAL/post_probs/fastml_sp/", gene_name + "_fastml_sp_char_indel.txt")
	else:
	    fastml_anc_seqs_marg_sp = join(fastml_species, "seq.marginal.txt") 
	    fastml_pp_sp = join("~/asr/FINAL/post_probs/fastml_sp/", gene_name + "_fastml_sp_char_only.txt")

	fastml_anc_seqs_joint_sp = join(fastml_species, "seq.joint.txt")

	# Read fastML trees with checks for existence
	if isfile(join(fastml_free, "tree.newick.txt")):
	    fastml_free_tree = Phylo.read(open(join(fastml_free, "tree.newick.txt")), "newick")
	else:
	    print(f"Warning: The file {join(fastml_free, 'tree.newick.txt')} does not exist.")

	if isfile(join(fastml_species, "tree.newick.txt")):
	    fastml_sp_tree = Phylo.read(open(join(fastml_species, "tree.newick.txt")), "newick")
	else:
	    print(f"Warning: The file {join(fastml_species, 'tree.newick.txt')} does not exist.")
	
	
	for x in fastml_free_tree.get_terminals() :
		if x.name == focal_name:
			scer_terminal = x
			break
	fastml_ancs_free = [fastml_free_tree.root]+[x for x in fastml_free_tree.root.get_path(scer_terminal)[:-1]]
	
	for x in fastml_sp_tree.get_terminals() :
		if x.name == focal_name:
			scer_terminal = x
			break
	fastml_ancs_sp = [fastml_sp_tree.root]+[x for x in fastml_sp_tree.root.get_path(scer_terminal)[:-1]]
	
	
	best_orf_list_anc = {}
	
	
	
	#######################################
	#######################################
	
	for x in prank_anc_tree_free.get_terminals() :
		if x.name == focal_name:
			scer_terminal = x
			break
	prank_ancs_free = [prank_anc_tree_free.root]+[x for x in prank_anc_tree_free.root.get_path(scer_terminal)[:-1]]
	
	for x in prank_anc_tree_sp.get_terminals() :
		if x.name == focal_name:
			scer_terminal = x
			break
	prank_ancs_sp = [prank_anc_tree_sp.root]+[x for x in prank_anc_tree_sp.root.get_path(scer_terminal)[:-1]]
	
	
	
	##################################
	###### fastml species joint ######
	##################################
	
	fastml_sp_joint_SEQS = [None for x in range(len(fastml_ancs_sp))]
	fastml_sp_just_names = [x.name for x in fastml_ancs_sp] 
	for Seq in SeqIO.parse(open(fastml_anc_seqs_joint_sp), 'fasta') :
		if Seq.id in fastml_sp_just_names :
			fastml_sp_joint_SEQS[fastml_sp_just_names.index(Seq.id)] = Seq	

	temp_orfs=[]
			
	for i in range(len(fastml_sp_joint_SEQS)) :
		anc = fastml_sp_joint_SEQS[i]
		anc.seq._data = ''.join(filter(lambda x : x!='-' , anc.seq._data))
		max_rand_orf = randomize_and_getcfs(scer_seq, anc, gene_name, no_reps)
		alignment = pairwise2.align.globalds(scer_seq, anc.seq._data, subs_mat, -3, -.1, one_alignment_only=True)
		res  = get_best_orf(alignment, scer_seq, anc, gene_name)
		temp_orfs.append(res[1])
		ranks = round(1-(np.array(max_rand_orf)<res[0]).sum() / len(max_rand_orf),4)
			
		terms = '_'.join(sorted(set([x.name for x in fastml_ancs_sp[i].get_terminals()])))
		
		print(gene_name, "FastML_joint", "species", fastml_sp_just_names[i], i, terms, round(res[0],3), ranks)		
		
	print(gene_name, "FastML_joint", "species", focal_name+"_spec", "NA", focal_name, -1, -1)
	best_orf_list_anc["FastML_joint_species"] = temp_orfs
	
	
	#####################################
	###### fastml species marginal ######
	#####################################
	
	fastml_sp_marg_SEQS = [None for x in range(len(fastml_ancs_sp))]
	fastml_sp_just_names = [x.name for x in fastml_ancs_sp] 
	for Seq in SeqIO.parse(open(fastml_anc_seqs_marg_sp), 'fasta') :
		if Seq.id in fastml_sp_just_names :
			fastml_sp_marg_SEQS[fastml_sp_just_names.index(Seq.id)] = Seq	

	temp_orfs=[]
			
	for i in range(len(fastml_sp_marg_SEQS)) :
		anc = fastml_sp_marg_SEQS[i]
		anc.seq._data = ''.join(filter(lambda x : x!='-' , anc.seq._data))
		max_rand_orf = randomize_and_getcfs(scer_seq, anc, gene_name, no_reps)
		alignment = pairwise2.align.globalds(scer_seq, anc.seq._data, subs_mat, -3, -.1, one_alignment_only=True)
		res  = get_best_orf(alignment, scer_seq, anc, gene_name)
		temp_orfs.append(res[1])
		ranks = round(1-(np.array(max_rand_orf)<res[0]).sum() / len(max_rand_orf),4)
			
		terms = '_'.join(sorted(set([x.name for x in fastml_ancs_sp[i].get_terminals()])))
		
		print(gene_name, "FastML_marginal", "species", fastml_sp_just_names[i], i, terms, round(res[0],3), ranks)		
		
	print(gene_name, "FastML_marginal", "species", focal_name+"_spec", "NA", focal_name, -1, -1)
	best_orf_list_anc["FastML_marginal_species"] = temp_orfs
	
	
	
	##################################
	###### fastml free joint ######
	##################################
	
	
	fastml_free_joint_SEQS = [None for x in range(len(fastml_ancs_free))]
	fastml_free_just_names = [x.name for x in fastml_ancs_free] 
	for Seq in SeqIO.parse(open(fastml_anc_seqs_joint_free), 'fasta') :
		if Seq.id in fastml_free_just_names :
			fastml_free_joint_SEQS[fastml_free_just_names.index(Seq.id)] = Seq	

	temp_orfs=[]
			
	for i in range(len(fastml_free_joint_SEQS)) :
		anc = fastml_free_joint_SEQS[i]
		anc.seq._data = ''.join(filter(lambda x : x!='-' , anc.seq._data))
		max_rand_orf = randomize_and_getcfs(scer_seq, anc, gene_name, no_reps)
		alignment = pairwise2.align.globalds(scer_seq, anc.seq._data, subs_mat, -3, -.1, one_alignment_only=True)
		res  = get_best_orf(alignment, scer_seq, anc, gene_name)
		temp_orfs.append(res[1])
		ranks = round(1-(np.array(max_rand_orf)<res[0]).sum() / len(max_rand_orf),4)
			
		terms = '_'.join(sorted(set([x.name for x in fastml_ancs_free[i].get_terminals()])))
		
		print(gene_name, "FastML_joint", "free", fastml_free_just_names[i], i, terms, round(res[0],3), ranks)		
		
	print(gene_name, "FastML_joint", "free", focal_name+"_spec", "NA", focal_name, -1, -1)
	best_orf_list_anc["FastML_joint_free"] = temp_orfs
	
	
	
	#####################################
	###### fastml free marginal ######
	#####################################
	
	fastml_free_marg_SEQS = [None for x in range(len(fastml_ancs_free))]
	fastml_free_just_names = [x.name for x in fastml_ancs_free] 
	for Seq in SeqIO.parse(open(fastml_anc_seqs_marg_free), 'fasta') :
		if Seq.id in fastml_free_just_names :
			fastml_free_marg_SEQS[fastml_free_just_names.index(Seq.id)] = Seq	

	temp_orfs=[]
			
	for i in range(len(fastml_free_marg_SEQS)) :
		anc = fastml_free_marg_SEQS[i]
		anc.seq._data = ''.join(filter(lambda x : x!='-' , anc.seq._data))
		max_rand_orf = randomize_and_getcfs(scer_seq, anc, gene_name, no_reps)
		alignment = pairwise2.align.globalds(scer_seq, anc.seq._data, subs_mat, -3, -.1, one_alignment_only=True)
		res  = get_best_orf(alignment, scer_seq, anc, gene_name)
		temp_orfs.append(res[1])
		ranks = round(1-(np.array(max_rand_orf)<res[0]).sum() / len(max_rand_orf),4)
			
		terms = '_'.join(sorted(set([x.name for x in fastml_ancs_free[i].get_terminals()])))
		
		print(gene_name, "FastML_marginal", "free", fastml_free_just_names[i], i, terms, round(res[0],3), ranks)		
		
	print(gene_name, "FastML_marginal", "free", focal_name+"_spec", "NA", focal_name, -1, -1)
	best_orf_list_anc["FastML_marginal_free"] = temp_orfs
	
	
	###########################
	###### prank species ######
	###########################
	
	prank_ancs_sp_SEQS = [None for x in range(len(prank_ancs_sp))]
	prank_ancs_sp_just_names = [x.name for x in prank_ancs_sp] 
	for Seq in SeqIO.parse(open(prank_anc_seqs_sp), 'fasta') :
		if Seq.id in prank_ancs_sp_just_names :
			prank_ancs_sp_SEQS[prank_ancs_sp_just_names.index(Seq.id)] = Seq	

	temp_orfs=[]
			
	for i in range(len(prank_ancs_sp_SEQS)) :
		anc = prank_ancs_sp_SEQS[i]
		anc.seq._data = ''.join(filter(lambda x : x!='-' , anc.seq._data))
		max_rand_orf = randomize_and_getcfs(scer_seq, anc, gene_name, no_reps)
		alignment = pairwise2.align.globalds(scer_seq, anc.seq._data, subs_mat, -3, -.1, one_alignment_only=True)
		res  = get_best_orf(alignment, scer_seq, anc, gene_name)
		temp_orfs.append(res[1])
		ranks = round(1-(np.array(max_rand_orf)<res[0]).sum() / len(max_rand_orf),4)
			
		terms = '_'.join(sorted(set([x.name for x in prank_ancs_sp[i].get_terminals()])))
		
		print(gene_name, "PRANK", "species", prank_ancs_sp_just_names[i], i, terms, round(res[0],3), ranks)		
		
	print(gene_name, "PRANK", "species", focal_name+"_spec", "NA", focal_name, -1, -1)
	best_orf_list_anc["PRANK_species"] = temp_orfs

	
	###########################
	###### prank free ######
	###########################
	
	prank_ancs_free_SEQS = [None for x in range(len(prank_ancs_free))]
	prank_ancs_free_just_names = [x.name for x in prank_ancs_free] 
	for Seq in SeqIO.parse(open(prank_anc_seqs_free), 'fasta') :
		if Seq.id in prank_ancs_free_just_names :
			prank_ancs_free_SEQS[prank_ancs_free_just_names.index(Seq.id)] = Seq	

	temp_orfs=[]
			
	for i in range(len(prank_ancs_free_SEQS)) :
		anc = prank_ancs_free_SEQS[i]
		anc.seq._data = ''.join(filter(lambda x : x!='-' , anc.seq._data))
		max_rand_orf = randomize_and_getcfs(scer_seq, anc, gene_name, no_reps)
		alignment = pairwise2.align.globalds(scer_seq, anc.seq._data, subs_mat, -3, -.1, one_alignment_only=True)
		res  = get_best_orf(alignment, scer_seq, anc, gene_name)
		temp_orfs.append(res[1])
		ranks = round(1-(np.array(max_rand_orf)<res[0]).sum() / len(max_rand_orf),4)
			
		terms = '_'.join(sorted(set([x.name for x in prank_ancs_free[i].get_terminals()])))
		
		print(gene_name, "PRANK", "free", prank_ancs_free_just_names[i], i, terms, round(res[0],3), ranks)		
		
	print(gene_name, "PRANK", "free", focal_name+"_spec", "NA", focal_name, -1, -1)
	best_orf_list_anc["PRANK_free"] = temp_orfs
	
	
	
	###########################
	##### prequel species #####
	###########################
	
			
	try:
		tree = Phylo.read(open(sp_tree_for_preq), "newick")
	except FileNotFoundError:
		print(f"Warning: Phylogenetic tree file '{sp_tree_for_preq}' not found. Skipping Prequel Species analysis.")
		return

	scer_terminal = None
	for x in tree.get_terminals() :
		if x.name == focal_name:
			scer_terminal = x
			break

	if scer_terminal is None:
		print("Terminal branch for focal species not found in the prequel species tree. Skipping...")
		return
	scer_anc = [tree.root]+[x for x in tree.root.get_path(scer_terminal)[:-1]]
	
	
	temp_orfs=[]
	for i in range(len(scer_anc)) :
		anc = scer_anc[i]
		try:
			for Seq in SeqIO.parse(open(join(prequel_species, "ORF_alignement."+anc.name+".fa")), 'fasta') :
				anc_seq = Seq
				break
		except FileNotFoundError:
			print(f"Warning: ORF alignment file for '{anc.name}' not found. Skipping this ancestor.")
			continue
		anc_seq.seq._data = ''.join(filter(lambda x : x!='-' , anc_seq.seq._data))
		alignment = pairwise2.align.globalds(scer_seq, anc_seq.seq._data, subs_mat, -3, -.1, one_alignment_only=True)
		res  = get_best_orf(alignment, scer_seq, anc_seq, gene_name)
		temp_orfs.append(res[1])
		max_rand_orf = randomize_and_getcfs(scer_seq, anc_seq, gene_name, no_reps)
		ranks=round(1-(np.array(max_rand_orf)<res[0]).sum() / len(max_rand_orf),4)
	
		terms = '_'.join(sorted(set([x.name for x in scer_anc[i].get_terminals()])))
		
		print(gene_name, "Prequel", "species", scer_anc[i].name, i, terms, round(res[0],3), ranks)		
	
	print(gene_name, "Prequel", "species", focal_name+"_spec", "NA", focal_name, -1, -1)		
	best_orf_list_anc["PREQUEL_species"] = temp_orfs	
	
	
	
	############################
	####### prequel free #######
	############################
	
		
	try:
		tree = Phylo.read(open(free_tree_for_preq), "newick")
	except FileNotFoundError:
		print(f"Warning: Phylogenetic tree file '{free_tree_for_preq}' not found. Skipping Prequel Free analysis.")
		return
	#tree = rename_branches(tree)
	scer_terminal = None
	for x in tree.get_terminals() :
		if x.name == focal_name:
			scer_terminal = x
			break
	if scer_terminal is None:
		print("Terminal branch for focal species not found in the prequel free tree. Skipping...")
		return

	scer_anc = [tree.root]+[x for x in tree.root.get_path(scer_terminal)[:-1]]
	
	temp_orfs=[]
	for i in range(len(scer_anc)) :
		anc = scer_anc[i]
		try:
			for Seq in SeqIO.parse(open(join(prequel_free, "ORF_alignement."+anc.name+".fa")), 'fasta') :
				anc_seq = Seq
				break
		except FileNotFoundError:
			print(f"Warning: ORF alignment file for '{anc.name}' not found. Skipping this ancestor.")
			continue
		anc_seq.seq._data = ''.join(filter(lambda x : x!='-' , anc_seq.seq._data))
		alignment = pairwise2.align.globalds(scer_seq, anc_seq.seq._data, subs_mat, -3, -.1, one_alignment_only=True)
		res  = get_best_orf(alignment, scer_seq, anc_seq, gene_name)
		temp_orfs.append(res[1])
		max_rand_orf = randomize_and_getcfs(scer_seq, anc_seq, gene_name, no_reps)
		ranks=round(1-(np.array(max_rand_orf)<res[0]).sum() / len(max_rand_orf),4)
	
		terms = '_'.join(sorted(set([x.name for x in scer_anc[i].get_terminals()])))
		
		print(gene_name, "Prequel", "free", scer_anc[i].name, i, terms, round(res[0],3), ranks)		
	
	print(gene_name, "Prequel", "free", focal_name+"_spec", "NA", focal_name, -1, -1)		
	best_orf_list_anc["PREQUEL_free"] = temp_orfs


	
	# Ensure the final output file is written
	orf_out = open(join(args[1], f"{args[1]}_bestRFC_ORFs.fasta"), 'w')
	for i in best_orf_list_anc.keys():
		for j in best_orf_list_anc[i] : 
			if type(j).__name__ == 'SeqRecord':
				j.id = i+"_"+j.id
				SeqIO.write(j, orf_out, "fasta")
	if comp_scer_acc:
	    SeqIO.write(comp_scer_acc, orf_out, "fasta")
	orf_out.close()	


if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
