# Final Group 12

## From Student Lead Tutorial 7: Genome Annotation Using Prokka, We cannot use Prokka because it is designed for prokaryotic genomes, instead we will use TransDecoder (predicting ORFs) and BLAST+ (annotation).




### Make sure you go to your ocean folder
``` bash
myocean
```
### Return to the Cloned Repository

``` bash
ls Student-Led-Tutorial-4
ls tutorial4
```

### Install and Load TransDecoder
```
module load anaconda3
conda create -n transdecoder_env -c bioconda transdecoder
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
module load cd-hit  # or conda install -c bioconda cd-hit
```

```
# Cluster within mock
cd-hit -i mock_trinity_out.Trinity.fasta.transdecoder.pep -o mock_cdhit.pep -c 0.95 -n 5

# Cluster within covid
cd-hit -i covid_trinity_out.Trinity.fasta.transdecoder.pep -o covid_cdhit.pep -c 0.95 -n 5

# Combine and find shared clusters
cat mock_cdhit.pep covid_cdhit.pep > combined.pep
cd-hit -i combined.pep -o combined_cdhit.pep -c 0.95 -n 5
```

### Using AlphaFold2 to predict 3D Structure (via ColabFold)

Useful Tutorial for AlphaFold2 and ChimeraX: https://www.youtube.com/watch?v=eLy7PdzRgLs 

Link for AlphaFold2:
https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb 





