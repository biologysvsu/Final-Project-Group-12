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
cd-hit -i mock_trinity_out.Trinity.fasta.transdecoder.pep -o mock_cdhit.pep -c 0.95 -n 5

# Cluster within covid
cd-hit -i covid_trinity_out.Trinity.fasta.transdecoder.pep -o covid_cdhit.pep -c 0.95 -n 5

# Combine and find shared clusters
cat mock_cdhit.pep covid_cdhit.pep > combined.pep
cd-hit -i combined.pep -o combined_cdhit.pep -c 0.95 -n 5
```

```
nano parse_cd_hit_clusters.py
```

### Paste and Run Python Script

```
from collections import defaultdict

clstr_file = "combined_cdhit.pep.clstr"
clusters = []
current_cluster = []

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

summary = defaultdict(list)

for idx, cluster in enumerate(clusters):
    covid_hits = 0
    mock_hits = 0
    seq_ids = []

    for line in cluster:
        parts = line.split()
        seq_id = parts[2].strip('>').strip('...')  # Extract sequence ID
        seq_ids.append(seq_id)

        if "covid" in seq_id.lower():
            covid_hits += 1
        elif "mock" in seq_id.lower():
            mock_hits += 1

    if covid_hits > 0 and mock_hits > 0:
        summary["shared"].append((idx, seq_ids))
    elif covid_hits > 0:
        summary["covid_only"].append((idx, seq_ids))
    elif mock_hits > 0:
        summary["mock_only"].append((idx, seq_ids))

# Save outputs
for category in ["shared", "covid_only", "mock_only"]:
    with open(f"{category}_ids.txt", "w") as out:
        for cluster_id, ids in summary[category]:
            for sid in ids:
                out.write(sid + "\n")
```
### Then run it
python parse_clusters.py

### Extract sequences with seqtk (install it first if using conda)
conda install -c bioconda seqtk

### Then run to extract sequences:

seqtk subseq combined_cdhit.pep shared_ids.txt > shared_seqs.fasta
seqtk subseq combined_cdhit.pep covid_only_ids.txt > covid_only_seqs.fasta
seqtk subseq combined_cdhit.pep mock_only_ids.txt > mock_only_seqs.fasta




### Using AlphaFold2 to predict 3D Structure (via ColabFold)

Useful Tutorial for AlphaFold2 and ChimeraX: https://www.youtube.com/watch?v=eLy7PdzRgLs 

Link for AlphaFold2:
https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb 





