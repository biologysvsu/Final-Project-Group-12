# Final-Project-Group-12
# Can Bioinformatics tools be used to visualize the structure of an unknown/hypothetical protein, and compare function of known functional proteins?
- SRA accession numbers:
- NP_659417, NP_689974, NP_001372653, NP_001372659, NP_001372657.
Software and Tools:
- BLAST, Phylogenetic analysis, Alpha fold
# Covid proteins were obtained from a denovo transcript assembly.
### File contains truncated and complete proteins. We will keep complete proteins only
```
mv covid_only_seqs.fasta covid.prot.fasta
```

```
awk '/^>/{keep=($0 ~ /complete/)} keep' covid.prot.fasta > complete_orfs.fa
```

# Lets blast to see whether we have any unknowns
1. Download the swissprot database:
```
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
```
2. Unzip it
```
gunzip uniprot_sprot.fasta.gz
mv uniprot_sprot.fasta swissprot.fasta
```
3. Make balst database
```
makeblastdb -in uniprot_sprot.fasta -dbtype prot -out swissprot
```
4. Blast
```
blastp -query complete_orfs.fa \
       -db swissprot \
       -evalue 1e-5 \
       -outfmt 6 \
       -max_target_seqs 1 \
       -num_threads 8 \
       -out complete_orfs_blast.tsv
```

# Lets extract the proteins that got not hit whatsoever, these maybe uncharacterized proteins.
1. extract no-hits
```
cut -f1 complete_orfs_blast.tsv | sort | uniq > with_hits.txt
grep "^>" complete_orfs.fa | cut -d' ' -f1 | sed 's/^>//' | sort > all_ids.txt
comm -23 all_ids.txt with_hits.txt > no_hit_ids.txt
seqtk subseq complete_orfs.fa no_hit_ids.txt > no_hit_orfs.fa
```
2. 
```
 awk '
  /^>/ {
    if (seq && length(seq) >= 120) {
      print header; print seq
    }
    header = $0
    seq = ""
    next
  }
  { seq = seq $0 }
  END {
    if (seq && length(seq) >= 120) {
      print header; print seq
    }
  }
' no_hit_orfs.fa > alphafold_nohit_candidates.fa
```
2. Count how many you got
```
grep -c ">" alphafold_nohit_candidates.fa
```
# Lets obtain hits with low percent identities under 30% similarity

```
awk '$3 < 30 {print $1 "\t" $3}' complete_orfs_blast.tsv > pid_map.tsv
cut -f1 pid_map.tsv | sort | uniq > low_pid_ids.txt
seqtk subseq complete_orfs.fa low_pid_ids.txt > tmp_low_pid.fa
```
- Change headers
```
awk '
BEGIN {
  FS = "\t"
  while ((getline < "pid_map.tsv") > 0) {
    pid[$1] = $2
  }
}
/^>/ {
  match($0, />[^ ]+/)
  id = substr($0, RSTART+1, RLENGTH-1)
  if (id in pid) {
    print ">" id " (identity: " pid[id] "%)"
  } else {
    print
  }
  next
}
{ print }
' tmp_low_pid.fa > alphafold_lowpid_candidates.fa
```
# Lets obtain proteins that hit proteins in the swissprot database that have been labeled uncharacterized
- Create joined_with_descriptions.tsv file
```
cut -f1,2 complete_orfs_blast.tsv | cut -f2 -d'|' --output-delimiter=' ' > blast_ids_and_hits.tsv
paste <(cut -f1 complete_orfs_blast.tsv) blast_ids_and_hits.tsv > blast_ids_and_hits_full.tsv

grep "^>" uniprot_sprot.fasta \
  | sed 's/^>//' \
  | awk -F'|' '{split($3,a," "); print $2, a[1], substr($0, index($0,a[1]))}' > swissprot_annotations.tsv

sort -k2 blast_ids_and_hits_full.tsv > sorted_hits.tsv
sort swissprot_annotations.tsv > sorted_annot.tsv

join -1 2 -2 1 sorted_hits.tsv sorted_annot.tsv > joined_with_descriptions.tsv
```

2. Get a list of hits
```
awk 'BEGIN{IGNORECASE=1} $0 ~ /uncharacterized|hypothetical|putative|virus|DUF/' joined_with_descriptions.tsv > weak_annotation_hits.tsv
```
3. Extract sequences
```
cut -d' ' -f2 weak_annotation_hits.tsv | sort | uniq > weak_annot_ids.txt

seqtk subseq complete_orfs.fa weak_annot_ids.txt > alphafold_weakannot_candidates.fa
```

# All files are in the repository
- alphafold_lowpid_candidates.fa has 2 sequences
- alphafold_nohit_candidates.fa has 77 sequences
- alphafold_weakannot_candidates.fa has 7 sequences
- alphafold_all_candidates.fa has has 86 sequences to work with

#Installing AlphaFold2 thru ColabFold

```
module load anaconda3
conda create -n colabfold -c conda-forge -c bioconda colabfold -y
conda activate colabfold
```
- Run alphafold, make sure your protein file is calles "alphafold_all_candidates.fa"
- Remove stop codons and fix name
```
awk '
/^>/ {
  match($0, /TRINITY_DN[0-9]+/); id = substr($0, RSTART, RLENGTH);
  match($0, /len:[0-9]+/); len = substr($0, RSTART+4, RLENGTH-4);
  print ">" id "_len" len;
  next
}
{
  gsub("\\*", ""); print
}
' alphafold_all_candidates.fa > alphafold_clean_named.fa

```
-create slurm script
```
vi colabfold.slurm
```
- type i and enter:
```
#!/bin/bash
#SBATCH --job-name=colabfold_run
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=8:00:00
#SBATCH --output=colabfold_%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jparedes@email.edu


# Load conda module
module load anaconda3

# Activate ColabFold environment
conda activate colabfold

# Run prediction
colabfold_batch alphafold_clean_named.fa out_dir
```
- run script
```
sbatch colabfold.slurm
```

# Install Folkseek to compare models to database and predict function.
```
module load anaconda3
conda create -n foldseek -c bioconda -c conda-forge foldseek
conda activate foldseek
```
