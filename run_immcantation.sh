AssignGenes.py igblast -s filtered_contig.fasta -b /usr/local/share/igblast --organism mouse --loci ig --format blast
MakeDb.py igblast -i filtered_contig_igblast.fmt7 -s filtered_contig.fasta -r /usr/local/share/germlines/imgt/mouse/vdj/imgt_mouse_IG*.fasta --10x filtered_contig_annotations.csv --extended
#MakeDb.py igblast -i filtered_contig_igblast.fmt7 -s filtered_contig.fasta -r /usr/local/share/igblast/fasta/imgt_mouse_ig_* --10x filtered_contig_annotations.csv --extended
ParseDb.py select -d filtered_contig_igblast_db-pass.tsv -f locus -u "IGH" --logic all --regex --outname heavy
ParseDb.py select -d filtered_contig_igblast_db-pass.tsv -f locus -u "IG[KL]" --logic all --regex --outname light
DefineClones.py -d heavy_parse-select.tsv --act set --model ham --norm len --dist 0.09
light_cluster.py -d heavy_parse-select_clone-pass.tsv -e light_parse-select.tsv -o filtered_contig_heavy_clone-light.tsv
CreateGermlines.py -d filtered_contig_heavy_clone-light.tsv -g dmask --cloned -r /usr/local/share/germlines/imgt/mouse/vdj/imgt_mouse_IGH*.fasta
BuildTrees.py -d filtered_contig_heavy_clone-light_germ-pass.tsv --minseq 2 --clean all --igphyml --nproc 12 --asr 0.9
