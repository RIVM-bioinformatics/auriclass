# Methodology

Under the hood, AuriClass conceptually works as follows:

``` mermaid
flowchart TD
  A[Argument parsing & input validation] -->|FastA data| B_fa[Create FastaAuriclass object];
  A --> |FastQ data| B_fq[Create FastqAuriclass object];
  B_fa --> C_fa[Sketch input FastA];
  B_fq --> C_fq[Sketch input FastQ];
  C_fa --> D[Run mash dist & select closest sample];
  C_fq --> D;
  D --> genome_size["Check whether genome size
                    is within expected range"];
  D --> E[Check if sample is close enough to any reference sample];
  E --> |Too distant| F_fail["FAIL: not Candida auris
                                Further QC skipped"];
  E --> |Close enough| F[Check is sample is closest to C. auris];
  F --> |Too distant| G_fail[WARN: other Candida/CUG-Ser1 clade sp.];
  F --> |Close enough| G[Select closest clade];
  G --> H[Final quality checks];
  G --> output
  H --> output

  genome_size --> output[Output]
  F_fail --> output
  G_fail --> output
```

!!! note
    Whether input files are FastQ or FastA is guessed by default based on `pyfastx` parsing. This can be overridden by specifying `--fastq` or `--fasta`.
    
    Steps that differ between FastQ and FastA data:
    
    - Sketching is performed differently: filtering for minimal kmer coverage is required for FastQ.
    - Genome size is estimated for FastQ using `mash sketch`, and parsed from FastA using `pyfastx`

`FastqAuriclass` and `FastaAuriclass` classes are based on `BasicAuriclass` which handles most of the attribute setting. The specific classes define some format-specific methods and how all methods should be called by the `run()` method. 
