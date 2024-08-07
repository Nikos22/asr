#!/bin/bash

# this script contains the commands used to run the phylogenetic reconstruction and ASR tools for the analyses of Vakirlis et al. GBE, 2024 .

INDIR=$1
SPECIES_TREE=$2
ORF_LIST=$3

cat $ORF_LIST | while read f ;
do
	SPECIES_TREE=$2
	#NAME="$(basename $f)" ;
	NAME=$f
	echo $NAME ;
	mkdir $NAME ;
	cd $NAME :
	#mafft $f > ${NAME}"_MAFFTali.fa" ;
	cp ${1}/${NAME}/MSA.MUSCLE.aln ${1}/${NAME}/MSA.MUSCLE_Seqs.Codes . ;
	cp MSA.MUSCLE.aln temptemp ;
	cat MSA.MUSCLE_Seqs.Codes | while read line ; do RN=$(echo $line | cut -f 1 -d ' ') ; NUM=$(echo $line | cut -f 2 -d ' ') ; cat temptemp | sed "s/${NUM}/${RN}/" > temptemptemp ; mv temptemptemp temptemp ; done ;
	mv temptemp ${NAME}"_ali.fa" ;
	grep ">"  ${NAME}"_ali.fa" | tr -d '>' > spNames.txt ;
	NO_OF_SP=$(wc -l spNames.txt | cut -f 1 -d ' ') ;
	echo $NO_OF_SP$SPECIES_TREE;
	if [ $NO_OF_SP -lt 8 ] ;
	then
		echo "No pruning..." ;
		gotree_amd64_linux prune -f spNames.txt -i ${SPECIES_TREE} -o pruned_species_tree.nwk -r ;
		SPECIES_TREE="pruned_species_tree.nwk"
	fi ;
	perl ~/asr/reconstructions_ORFs/Fasta2Phylip.pl ${NAME}"_ali.fa" ${NAME}"_ali.phy" ;
	
	## estimate trees, using GTR (REV) model, 4 gamma rate cats and empirical base freqs
	raxml-ng --seed 12546582 --evaluate --msa ${NAME}"_ali.fa" --model GTR+F+G --threads 2 --tree $SPECIES_TREE --prefix ${NAME}"_fixedTop" ;
	raxml-ng --seed 12546582 --msa ${NAME}"_ali.fa"  --model GTR+F+G --prefix ${NAME}"_freeTop" --threads 2 ;
	
	~/bin/gotree_amd64_linux reroot outgroup -i ${NAME}"_fixedTop.raxml.bestTree" Seub Suva > ${NAME}"_fixedTop.raxml.bestTree_CORPAR" ;
	~/bin/gotree_amd64_linux reroot midpoint -i ${NAME}"_freeTop.raxml.bestTree" > ${NAME}"_freeTop.raxml.bestTree_midRoot" ;
	
	
	phyloFit --seed 12546582 --tree ${NAME}"_freeTop.raxml.bestTree_midRoot" --msa ${NAME}"_ali.fa"  --out-root ${NAME}"_freeTopology.phylofit"
	phyloFit --seed 12546582 --tree ${NAME}"_fixedTop.raxml.bestTree_CORPAR" --msa ${NAME}"_ali.fa"  --out-root ${NAME}"_speciesTopology.phylofit"
	
	
	grep TREE ${NAME}"_freeTopology.phylofit.mod" | sed 's/^TREE\: //' > temp.tree ;
	~/bin/gotree_amd64_linux reroot midpoint -i temp.tree | ~/bin/gotree_amd64_linux rename -l 1 --tips=false --internal -a | sed 's/N000/N/g' > temp.tree.midroot ;	
	
	TR=$(cat temp.tree.midroot) ;
	cat  ${NAME}"_freeTopology.phylofit.mod" | sed "s/^TREE\:.*/TREE: $TR/" > ${NAME}"_freeTopology.phylofit_corTree.mod" ;
	
	grep TREE ${NAME}"_speciesTopology.phylofit.mod" | sed 's/^TREE\: //' > temp.tree ;
	~/bin/gotree_amd64_linux reroot outgroup -i temp.tree Seub Suva | ~/bin/gotree_amd64_linux rename -l 1 --tips=false --internal -a | sed 's/N000/N/g' > temp.tree.corout ;

	TR=$(cat temp.tree.corout) ;
        cat  ${NAME}"_speciesTopology.phylofit.mod" | sed "s/^TREE\:.*/TREE: $TR/" > ${NAME}"_speciesTopology.phylofit_corTree.mod" ;
	
	mkdir prank_freeTop prank_spTop prequel_freeTop prequel_spTop fastml_freeTop fastml_spTop
	prank -d=${NAME}"_ali.fa"  -support -showall -keep -F -once -o=./prank_freeTop/ORF_alignment -t=${NAME}"_freeTop.raxml.bestTree_midRoot"
	prank -d=${NAME}"_ali.fa"  -support -showall -keep -F -once -o=./prank_spTop/ORF_alignment -t=${NAME}"_fixedTop.raxml.bestTree_CORPAR"
	
	cp ${NAME}"_ali.fa" fastml_freeTop ;
	cp ${NAME}"_ali.fa" fastml_spTop ;
	perl ~/asr/programs/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_File ${NAME}"_ali.fa" --seqType nuc --Tree ${NAME}"_freeTop.raxml.bestTree_midRoot" --SubMatrix GTR --OptimizeBL no --indelReconstruction ML --outDir ${PWD}/fastml_freeTop/
	perl ~/asr/programs/FastML.v3.11/www/fastml/FastML_Wrapper.pl --MSA_File ${NAME}"_ali.fa" --seqType nuc --Tree ${NAME}"_fixedTop.raxml.bestTree_CORPAR" --SubMatrix GTR --OptimizeBL no --indelReconstruction ML --outDir ${PWD}/fastml_spTop/
	
	prequel ${NAME}"_ali.fa" ${NAME}"_freeTopology.phylofit_corTree.mod" ./prequel_freeTop/ORF_alignement -n
	prequel ${NAME}"_ali.fa" ${NAME}"_speciesTopology.phylofit_corTree.mod" ./prequel_spTop/ORF_alignement -n
	prequel ${NAME}"_ali.fa" ${NAME}"_freeTopology.phylofit_corTree.mod" ./prequel_freeTop/ORF_alignement
        prequel ${NAME}"_ali.fa" ${NAME}"_speciesTopology.phylofit_corTree.mod" ./prequel_spTop/ORF_alignement
	cd ../ :
	
done
