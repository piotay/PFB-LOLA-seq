# PFB-RNA_seq-2024


# Github Setup

    - Create a fork of the https://github.com/rdalipo1/PFB-LOLA-seq repository


    - git clone the ssh for you fork


    - Make sure your fork is updated before pushing from your local


    - To contribute to the main, you need to push to your fork first, then contribute from the github website

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
