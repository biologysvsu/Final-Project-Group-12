# Final-Project-Group-12
# Can Bioinformatics tools be used to visualize the structure of an unknown/hypothetical protein, and compare function of known functional proteins?
- SRA accession numbers:
- NP_659417, NP_689974, NP_001372653, NP_001372659, NP_001372657.
Software and Tools:
- BLAST, Phylogenetic analysis, Alpha fold
# Covid proteins were obtained from a denovo transcript assembly.
### File contains truncated and complete proteins. We will keep complete proteins only
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
