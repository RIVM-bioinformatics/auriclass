# Input & output

## Input

AuriClass takes an arbitrary number of Fastq or Fasta files as input (however not Fastq and Fasta in a single analysis). AuriClass assumes all files originate from a single sample.

## Reference data (included)

Two files with reference data are needed to start the analysis. Default files are included:

| File | Description |
|---|---|
| `data/clade_config.csv` | CSV table listing two columns: filename and associated clade. All included species that are not *Candida auris* are listed as "outgroup". |
| `data/Candida_auris_clade_references.msh` | Mash sketch of reference genomes of *Candida auris* and related species. |

See [Reference genomes](reference_genomes.md) for more information on the default database and see the [FAQ](faq.md) for how to build your own database or change calibrated settings.

## Output

AuriClass outputs a tab-separated text file containing the main output, and a log file which save info written to STDERR.

The columns in the output report have the following information:

| Column | Description |
|---|---|
| Sample | Sample name. Default value is "isolate". |
| Clade | Predicted clade. If the closest sample is not *Candida auris*, this field will contain "not Candida auris". |
| Mash_distance_from_closest_reference | Mash distance from the closest reference genome. |
| QC_decision | "FAIL" is QC_species is "FAIL". "WARN" if any other QC step indicates "WARN". Otherwise "PASS". |
| QC_species | "FAIL" if the mash distance from the closest reference is more than `--non_candida_threshold`. This indicates that the sample is not closely related to any *Candida* spp. |
| QC_other_Candida | "WARN" if the closest reference genome is not a *Candida auris* reference genome. This indicates that the sample is closer to another species than *Candida auris*. |
| QC_genome_size | "WARN" if the estimated genome size is outside the specified range. This is probably not accurate in case you run a diploid sample. |
| QC_multiple_hits | "WARN" if at least one sample is within the 99% Mash error bounds of the closest hit (from `mash bounds -k [KMERSIZE] -p 0.99`). This should not happen if the appropriate settings and database are used. |
| QC_possible_new_clade | "WARN" if the closest reference genome is *Candida auris*, but the Mash distance is higher than the `--new_clade_threshold`. |