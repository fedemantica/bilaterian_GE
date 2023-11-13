wget http://cisbp.ccbr.utoronto.ca/data/2.00/DataFiles/Bulk_downloads/EntireDataset/TF_Information_all_motifs.txt.zip
wget http://cisbp.ccbr.utoronto.ca/data/2.00/DataFiles/Bulk_downloads/EntireDataset/PWMs.zip

unzip -p TF_Information_all_motifs.txt.zip | head -n 1  > TF_INFO_BILATERIA.txt
unzip -p TF_Information_all_motifs.txt.zip | grep -P "Homo_sapiens|Mus_musculus|Drosophila_melanogaster|Rattus_norvegicus|Danio_rerio|Tetraodon_nigroviridis|Xenopus_tropicalis|Gallus_gallus|Xenopus_laevis|Meleagris_gallopavo|Monodelphis_domestica|Strongylocentrotus_purpuratus|Bombyx_mori|Schistosoma_mansoni|Oikopleura_dioica|Oryzias_latipes|Anolis_carolinensis|Takifugu_rubripes|Ornithorhynchus_anatinus|Gasterosteus_aculeatus|Sus_scrofa|Anopheles_gambiae|Branchiostoma_floridae|Taeniopygia_guttata|Aedes_aegypti|Oryctolagus_cuniculus|Cavia_porcellus|Apis_mellifera|Oncorhynchus_tshawytscha|Choloepus_hoffmanni|Dasypus_novemcinctus|Ochotona_princeps|Cupiennius_salei|Daphnia_pulex|Pan_paniscus|Pan_troglodytes|Heterodontus_francisci|Hydra_magnipapillata|Pediculus_humanus|Petromyzon_marinus|Tursiops_truncatus|Tupaia_belangeri|Acyrthosiphon_pisum|Acipenser_baerii|Ailuropoda_melanoleuca|Trionyx_sinensis|Tarsius_syrichta|Canis_familiaris|Pongo_abelii|Procavia_capensis|Pteropus_vampyrus|Drosophila_yakuba|Drosophila_virilis|Drosophila_pseudoobscura|Microcebus_murinus|Marmota_monax"  | awk '{if($9=="D"){print $0;}}' OFS="\t" >> TF_INFO_BILATERIA.txt

unzip -d matrices PWMs.zip

mkdir -p perfam
tail -n +2 TF_INFO_BILATERIA.txt | perl -ne '{@vec=split("\t",$_); $ac=$vec[3];
   $id=$vec[6].":".$vec[9].":".$vec[10]; $cc=$vec[7];
   $family=$vec[9];
   $family =~ s/\//-/g;
   $family =~ s/\s/_/g;
   $family =~ s/,/_plus_/g;
   $file="matrices/pwms/".$ac.".txt";
   $fileOut="perfam/pwms.".$family.".tf";
   $n=0;
   open(T,$file) or die "File $file does not exist\n"; while(<T>){$n++;}close(T);
   if($n>1){
					open(OUT,">>",$fileOut) or die "Can not open $fileOut\n";
					print OUT "AC\t$ac\nXX\nID\t$id\nXX\nDE\t$cc\n";
					open(T,$file);while($l=<T>){print OUT $l;}close(T);
					print OUT "XX\nCC\tSpecies:$cc\nXX\n\/\/\n";
					close(OUT);}
   }'

sed -i 's/Pos/P0/g' perfam/*


#####
##
## For each family = ${FAM}
##
######

rsat convert-matrix -v 1 -i perfam/pwms.${FAM}.tf -from tf -to cb -o perfam/pwms.${FAM}.cb -decimals 6 -return counts
awk '{if($1~/>/)print $1;else print $0;}' perfam/pwms.${FAM}.cb > perfam/pwms.${FAM}.corrected.cb
gimme cluster -t 0.9999 perfam/pwms.${FAM}.corrected.cb gimme_cluster/${FAM}

gimme cluster -t 0.9999 -N 1 perfam/pwms_${FAM}.cb ${FAM}
for i in $(ls perfam/*.tf | sed 's/perfam\/pwms.//g' | sed 's/.tf//g')
do
	sed "s/>/>${i}_/g" gimme_cluster/${i}/clustered_motifs.pfm >> clustered_motifs.pfm #-i
	awk -v tf=$i '{  $1=tf"_"$1; print $0;}' OFS="\t" gimme_cluster/${i}/cluster_key.txt >> clustered_key_2Motif.txt
done

################## END ##################

for i in $(ls perfam/*.tf | sed 's/perfam\/pwms.//g' | sed 's/.tf//g')
do
	sed "s/>/>${i}_/g" gimme_cluster/${i}/clustered_motifs.pfm >> clustered_motifs.pfm #-i
	awk -v tf=$i '{  $1=tf"_"$1; print $0;}' OFS="\t" gimme_cluster/${i}/cluster_key.txt >> clustered_key_2Motif.txt
done


rsat convert-matrix convert-matrix  -from cluster-buster -to cluster-buster \
-i clustered_motifs.pfm 1 -multiply 1 -decimals 1 -perm 0 -bg_pseudo 0.01 \
-return counts -to cluster-buster | grep '>' | cut -f 1,6 -d ' ' | sed 's/>//g' | sed 's/\/size=//g' > clustered_motifs_sizes.txt
awk '{if ($2 < 5) print $1}' clustered_motifs_sizes.txt > clusters2remove.txt

perl -e 'open(LI,$ARGV[0]);while(<LI>){chomp $_;$h{$_}=1;}close(LI);
          $print=1;open(FA,$ARGV[1]);
          while(<FA>){
             if($_ =~ /^>(.+)/){ $l=$1;
                if($h{$l}==1){$print=0;}else{$print=1;}
             }
             if($print){print $_;
             }
          }close(FA);' clusters2remove.txt clustered_motifs.pfm > clustered_motifs_filtered_size4.pfm
          




