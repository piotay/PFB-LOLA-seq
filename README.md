![Alt text](https://raw.githubusercontent.com/piotay/PFB-LOLA-seq/refs/heads/main/pfb_lola_seq/scrna-seq/atac/LOLA-SEQ-10-28-2024.jpg)


# Github Setup

Create a fork of the https://github.com/rdalipo1/PFB-LOLA-seq repository
    
- Clone the ssh for you fork
```bash
git clone git@github.com:usr/fork_name.git
```
    
Make sure your fork is updated before pushing from your local	
- Sync on the website
- Pull on local
```git pull```

To contribute to the main, you need to push to your fork first, then contribute from the github website
- Push to your fork
```bash
git add
git commit -m 'message'
git push
```
- Contribute from fork github website

# Development Setup

Create and activate a Python virtual environment with mamba:
```bash
$ mamba create -n lola-seq python=3.12
...
$ mamba activate lola-seq
```

Move to the top-level project directory and perform an "editable" install of the project:
```bash
$ cd PFB-LOLA-seq
PFB-LOLA-seq $ pip install -e .
```

This will install the project dependencies in the `lola-seq` virtual environment.

# SRAToolkit Instructions

Download from mamba:
```bash
mamba activate lola-seq
mamba install -c bioconda -c conda-forge sra-tools=3.1.1
```

Instructions on how to use SRAtoolkit:
- Prefetch to grab the SRA file from NCBI: https://erilu.github.io/python-fastq-downloader/
- fastq-dump to convert SRA file to fastq: https://edwards.flinders.edu.au/fastq-dump/
