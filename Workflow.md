# Final Group 12

## From Student Lead Tutorial 7: Genome Annotation Using Prokka, We cannot use Prokka because it is designed for prokaryotic genomes, instead we will use TransDecoder (predicting ORFs) and BLAST+ (annotation).


## This is a continuation of Student-Lead Tutorial 4




### Make sure you go to your ocean folder
``` bash
myocean
```
### Return to the Cloned Repository

``` bash
cd Student-Led-Tutorial-4
cd tutorial4
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


### Swissprot database
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz


### Load blast
conda install -c bioconda blast


### Using AlphaFold2 to predict 3D Structure (via ColabFold)

Useful Tutorial for AlphaFold2 and ChimeraX: https://www.youtube.com/watch?v=eLy7PdzRgLs 

Link for AlphaFold2:
https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb 





