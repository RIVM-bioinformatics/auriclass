# Installation

## Conda

The easiest way to install AuriClass is through mamba/conda:

```
mamba create -n env_auriclass -c bioconda -c conda-forge auriclass
conda activate env_auriclass
```

## Containers

Alternatively, [biocontainers](https://quay.io/repository/biocontainers/auriclass) has built a container from the AuriClass bioconda package. This container can be used with different container software, e.g. Singularity:

```
singularity pull docker://quay.io/biocontainers/auriclass:0.5.3--pyhdfd78af_0
```

!!! note
    `auriclass:latest` is not defined in biocontainers, so make sure to specify the latest bioconda version and build.