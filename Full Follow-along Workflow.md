# Final-Project-Group-12
# Can Bioinformatics tools be used to visualize the structure of an unknown/hypothetical protein, and compare function of known functional proteins?
## **Dataset Description**
We will download RNA-seq data for mock-infected and COVID-19-infected human cell lines from the provided SRA project.
- **Mock-infected cells**:
  - SRA Accessions: `SRR11412215`, `SRR11412216`
- **COVID-19-infected cells**:
  - SRA Accessions: `SRR11412227`, `SRR11412228`
Software and Tools:
- SRAToolkit, FastQC, Trinity, BLAST, Phylogenetic analysis, Alpha fold
### Make sure you go to your ocean folder
``` bash
myocean
```
### Clone the Repository
``` bash
git clone
```
``` bash
cd Final-Project-Group-12
mkdir tutorial
cd tutorial
```
### **Download Instructions**
Use the SRA Toolkit (must be installed beforehand if run locally, otherwise available in most HPCs) to download paired-end FASTQ files:

   ``` bash
   # Example for a mock-infected sample. More replicates are always better, so repeat step for each SRA    accession.
# Example for a COVID-19-infected sample. More replicates are always better, so repeat step for each SRA accession.


module load sra-toolkit
prefetch --max-size 100G SRR11412215 #Handles large files efficiently (downloads in chunks to avoid corruption) 
fastq-dump --gzip SRR11412215

prefetch --max-size 100G SRR11412216
fastq-dump --gzip SRR11412216

prefetch --max-size 100G SRR11412227
fastq-dump --gzip SRR11412227

prefetch --max-size 100G SRR11412228
fastq-dump --gzip SRR11412228
```
### **Tasks and Deliverables**
#### **Part 1: Data Preparation**
1. Unzip the .fastq.qz files
   ```bash
   gunzip SRR11412215.fastq.gz
   gunzip SRR11412216.fastq.gz
   gunzip SRR11412227.fastq.gz
   gunzip SRR11412228.fastq.gz

1. Verify the quality of RNA-seq reads using FastQC (optional but recommended):
   ```bash
   module load FastQC
   mkdir fastqc_data_out
   fastqc -t 4 *.fastq -o fastqc_data_out
2. Run multiqc
```
module load python/3.8.6
multiqc --dirs fastqc_data_out --filename multiqc_raw_data.html
```
-Visualize multiqc file with QIIME 2.
```
git add multiqc_raw_data.html
git commit -m "multiqc data"
git push
```

3. Combine FASTQ files for mock-infected and COVID-19-infected groups:
   ```bash
   cat SRR11412215.fastq SRR11412216.fastq > mock_combined.fastq
   cat SRR11412227.fastq SRR11412228.fastq > covid_combined.fastq
   
#### **Part 2: Run Trinity (MEMORY DEMANDING!)**
1. Run Trinity for de novo transcriptome assembly:
```
salloc --mem=128G --cpus-per-task=128 --time=06:00:00

module load Trinity/2.15.1

Trinity --seqType fq --max_memory 128G \
        --single mock_combined.fastq \
        --CPU 128 --output mock_trinity_out
```

  - Replace mock_combined.fastq with covid_combined.fastq
```
Trinity --seqType fq --max_memory 128G \
        --single covid_combined.fastq \
        --CPU 128 --output covid_trinity_out
```

2. Output:
- Trinity.fasta: Contains assembled transcripts.


### Make sure you go to your ocean folder
``` bash
myocean
```
### Return to the Cloned Repository

``` bash
cd Final-Project-Group-12
cd tutorial
```

### Install and Load TransDecoder
```
module load anaconda3
conda create -n transdecoder_env -c bioconda -c conda-forge transdecoder cd-hit -y
conda activate transdecoder_env
```

### Running TransDecoder to Predict Coding Regions

Mock
```
TransDecoder.LongOrfs -t mock_trinity_out.Trinity.fasta
TransDecoder.Predict -t mock_trinity_out.Trinity.fasta
```
COVID
```
TransDecoder.LongOrfs -t covid_trinity_out.Trinity.fasta
TransDecoder.Predict -t covid_trinity_out.Trinity.fasta
```

### Using CD-HIT to find Unique & Shared Transcripts

```
# Cluster within mock
cd-hit -i combined.pep -o combined_cdhit.pep -c 0.95 -n 5 -d 0
```
```
cd-hit -i mock_trinity_out.Trinity.fasta.transdecoder.pep -o mock_cdhit.pep -c 0.95 -n 5 -d 0
cd-hit -i covid_trinity_out.Trinity.fasta.transdecoder.pep -o covid_cdhit.pep -c 0.95 -n 5 -d 0
cat mock_cdhit.pep covid_cdhit.pep > combined.pep
cd-hit -i combined.pep -o combined_cdhit.pep -c 0.95 -n 5 -d 0
```
```
# Include .p1/.p2 suffix in the IDs
grep "^>" mock_cdhit.pep | cut -d ' ' -f1 | sed 's/^>//' > mock_ids_raw.txt
grep "^>" covid_cdhit.pep | cut -d ' ' -f1 | sed 's/^>//' > covid_ids_raw.txt
```

```
grep -oP '(?<=\>).+?(?=\.\.\.)' combined_cdhit.pep.clstr | head
```
```
nano parse_clusters.py
```

### Paste and Run Python Script

```
import re
from collections import defaultdict

# Load mock and covid IDs
def load_ids(filename):
    with open(filename) as f:
        return set(line.strip() for line in f if line.strip())

mock_ids = load_ids("mock_ids_raw.txt")
covid_ids = load_ids("covid_ids_raw.txt")

clstr_file = "combined_cdhit.pep.clstr"
clusters = []
current_cluster = []

# Read clusters
with open(clstr_file, 'r') as f:
    for line in f:
        if line.startswith(">Cluster"):
            if current_cluster:
                clusters.append(current_cluster)
                current_cluster = []
        else:
            current_cluster.append(line.strip())
    if current_cluster:
        clusters.append(current_cluster)

# Analyze clusters
summary = defaultdict(list)

for idx, cluster in enumerate(clusters):
    cluster_ids = []
    has_mock = False
    has_covid = False

    for line in cluster:
        # Extract FASTA ID between '>' and '...'
        match = re.search(r'>(.+?)\.\.\.', line)
        if not match:
            continue
        seq_id = match.group(1)
        cluster_ids.append(seq_id)

        if seq_id in mock_ids:
            has_mock = True
        if seq_id in covid_ids:
            has_covid = True

    if has_mock and has_covid:
        summary["shared"].append((idx, cluster_ids))
    elif has_mock:
        summary["mock_only"].append((idx, cluster_ids))
    elif has_covid:
        summary["covid_only"].append((idx, cluster_ids))

# Write outputs
for category in ["shared", "mock_only", "covid_only"]:
    with open(f"{category}_ids.txt", "w") as out:
        written = set()
        for _, ids in summary[category]:
            for sid in ids:
                if sid not in written:
                    out.write(sid + "\n")
                    written.add(sid)

print("✅ Done! Wrote:")
print(" - shared_ids.txt")
print(" - mock_only_ids.txt")
print(" - covid_only_ids.txt")

```
### Then run it
```
python parse_clusters.py
```

### Check output 
```
wc -l mock_only_ids.txt covid_only_ids.txt shared_ids.txt
```
### Extract sequences with seqtk (install it first if using conda)
```
conda install -c bioconda seqtk
```
### Then run to extract sequences:
```
seqtk subseq combined_cdhit.pep shared_ids.txt > shared_seqs.fasta
seqtk subseq combined_cdhit.pep covid_only_ids.txt > covid_only_seqs.fasta
seqtk subseq combined_cdhit.pep mock_only_ids.txt > mock_only_seqs.fasta
```

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






