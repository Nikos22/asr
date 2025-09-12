#!/bin/bash

# this script contains the commands used to run the phylogenetic reconstruction and ASR tools for the analyses of Vakirlis et al. GBE, 2024 .

INPUT_SEQS=$1
SPECIES_TREE=$2
NO_OF_SPECIES_IN_TREE=$3
OUTGROUP_SPECIES="Seub Suva" # specify the tip names of the outgroup in your tree, separated by spaces 

NAME=${INPUT_SEQS%%\.*}
echo $NAME ;
mkdir $NAME ;
cp $INPUT_SEQS $SPECIES_TREE $NAME ;
cd $NAME ;

## align the sequences
mafft $INPUT_SEQS > ${NAME}"_ali.fa" ;

## get sequence names, count them and check if their number matches the tree
grep ">"  ${NAME}"_ali.fa" | tr -d '>' > spNames.txt ;
NO_OF_SP=$(wc -l spNames.txt | cut -f 1 -d ' ') ;

if [ $NO_OF_SP -lt $NO_OF_SPECIES_IN_TREE ] ;
then
	# if there are less sequences then prune the tree
	echo "Pruning species tree" ;
	gotree_amd64_linux prune -f spNames.txt -i ${SPECIES_TREE} -o pruned_species_tree.nwk -r ;
	SPECIES_TREE="pruned_species_tree.nwk"
fi ;

# run a custom script to convert to PHYLIP, change path accordingly
perl ../Fasta2Phylip.pl ${NAME}"_ali.fa" ${NAME}"_ali.phy" ;

## estimate trees, using GTR (REV) model, 4 gamma rate cats and empirical base freqs
raxml-ng --seed 12546582 --evaluate --msa ${NAME}"_ali.fa" --model GTR+F+G --threads 2 --tree $SPECIES_TREE --prefix ${NAME}"_fixedTop" ;
raxml-ng --seed 12546582 --msa ${NAME}"_ali.fa"  --model GTR+F+G --prefix ${NAME}"_freeTop" --threads 2 ;

# necessary adjustment to reroot the tree correctly due to removal of parentheses by RAXML, here is where the outgroups are needed
gotree_amd64_linux reroot outgroup -i ${NAME}"_fixedTop.raxml.bestTree" $OUTGROUP_SPECIES > ${NAME}"_fixedTop.raxml.bestTree_CORPAR" ;

# reroot at midpoint when going with a free topology
gotree_amd64_linux reroot midpoint -i ${NAME}"_freeTop.raxml.bestTree" > ${NAME}"_freeTop.raxml.bestTree_midRoot" ;

# tree preparation for PREQUEL
phyloFit --seed 12546582 --tree ${NAME}"_freeTop.raxml.bestTree_midRoot" --msa ${NAME}"_ali.fa"  --out-root ${NAME}"_freeTopology.phylofit"
phyloFit --seed 12546582 --tree ${NAME}"_fixedTop.raxml.bestTree_CORPAR" --msa ${NAME}"_ali.fa"  --out-root ${NAME}"_speciesTopology.phylofit"


grep TREE ${NAME}"_freeTopology.phylofit.mod" | sed 's/^TREE\: //' > temp.tree ;
gotree_amd64_linux reroot midpoint -i temp.tree | gotree_amd64_linux rename -l 1 --tips=false --internal -a | sed 's/N000/N/g' > temp.tree.midroot ;	

TR=$(cat temp.tree.midroot) ;
cat  ${NAME}"_freeTopology.phylofit.mod" | sed "s/^TREE\:.*/TREE: $TR/" > ${NAME}"_freeTopology.phylofit_corTree.mod" ;

grep TREE ${NAME}"_speciesTopology.phylofit.mod" | sed 's/^TREE\: //' > temp.tree ;
gotree_amd64_linux reroot outgroup -i temp.tree $OUTGROUP_SPECIES | gotree_amd64_linux rename -l 1 --tips=false --internal -a | sed 's/N000/N/g' > temp.tree.corout ;

TR=$(cat temp.tree.corout) ;
cat  ${NAME}"_speciesTopology.phylofit.mod" | sed "s/^TREE\:.*/TREE: $TR/" > ${NAME}"_speciesTopology.phylofit_corTree.mod" ;
	
# run PRANK
mkdir prank_freeTop prank_spTop prequel_freeTop prequel_spTop fastml_freeTop fastml_spTop
prank -d=${NAME}"_ali.fa"  -support -showall -keep -F -once -o=./prank_freeTop/ORF_alignment -t=${NAME}"_freeTop.raxml.bestTree_midRoot"
prank -d=${NAME}"_ali.fa"  -support -showall -keep -F -once -o=./prank_spTop/ORF_alignment -t=${NAME}"_fixedTop.raxml.bestTree_CORPAR"

## run FastML
cp ${NAME}"_ali.fa" fastml_freeTop ;
cp ${NAME}"_ali.fa" fastml_spTop ;

## change the path to FastML accordingly
perl ~/asr/programs/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_File ${NAME}"_ali.fa" --seqType nuc --Tree ${NAME}"_freeTop.raxml.bestTree_midRoot" --SubMatrix GTR --OptimizeBL no --indelReconstruction ML --outDir ${PWD}/fastml_freeTop/
perl ~/asr/programs/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_File ${NAME}"_ali.fa" --seqType nuc --Tree ${NAME}"_fixedTop.raxml.bestTree_CORPAR" --SubMatrix GTR --OptimizeBL no --indelReconstruction ML --outDir ${PWD}/fastml_spTop/
	
## run PREQUEL
prequel ${NAME}"_ali.fa" ${NAME}"_freeTopology.phylofit_corTree.mod" ./prequel_freeTop/ORF_alignement -n
prequel ${NAME}"_ali.fa" ${NAME}"_speciesTopology.phylofit_corTree.mod" ./prequel_spTop/ORF_alignement -n
prequel ${NAME}"_ali.fa" ${NAME}"_freeTopology.phylofit_corTree.mod" ./prequel_freeTop/ORF_alignement
prequel ${NAME}"_ali.fa" ${NAME}"_speciesTopology.phylofit_corTree.mod" ./prequel_spTop/ORF_alignement

cd ../ ;
