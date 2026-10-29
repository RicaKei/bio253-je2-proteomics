# data/ – local only, never committed

The data are **unpublished**. Everything in this folder except this README is
git-ignored. Do not remove those lines from `.gitignore`, and never upload the files
anywhere (GitHub, cloud storage, chatbots).

## How to get the data

1. Get the course data folder `SAJE2/` from the lecturer (course share).
2. From the repository root, run:

   ```bash
   Rscript scripts/00_fetch_data.R /path/to/SAJE2          # copies (~135 MB)
   Rscript scripts/00_fetch_data.R /path/to/SAJE2 --link   # symlinks instead of copies
   ```

   You can also set the environment variable `BIO253_DATA_DIR=/path/to/SAJE2`.

## What you get

```
data/raw/
  MANIFEST.tsv                       source path, size and md5 of every copied file
  total/DE_total.xlsx                prolfquapp DEA, total proteome (Oct 2026 run)
  phospho_STYH/DE_phospho_STYH.xlsx  prolfquapp DEA, phospho-sites, S/T/Y/H search (Oct 2026)
  phospho_STY/DE_phospho_STY.xlsx    prolfquapp DEA, phospho-sites, S/T/Y search (Aug 2026, older model)
  sites/combined_site_STY.tsv        FragPipe site tables (raw, wide)
  sites/combined_site_STYH.tsv
  resources/JE2_proteins.fasta       search database (JE2 CP020619.1 + pseudogene translations + contaminants)
  resources/JE2_genome.gff           genome annotation
  resources/string_*.txt             STRING v12 aliases, enrichment terms, physical links for JE2
data/derived/                        tidy long tables (tab-separated, gzipped)
  protein_annotation.tsv.gz          locus_tag, product, protein_id, location, NCTC8325 ortholog
  samples.tsv                        27 samples: genotype, medium, replicate, group
  dea_<dataset>.tsv.gz               DEA results, one row per protein/site x contrast
  abundances_<dataset>.tsv.gz        normalised abundances, one row per protein/site x sample
  sites_STY.tsv.gz, sites_STYH.tsv.gz  FragPipe sites, one row per site x sample (0 -> NA)
data/reference/                      lecturer reference tables, if available (e.g. pHis tiers)
```

`<dataset>` is `total`, `phospho_STYH` or `phospho_STY`. Load the tables with the
helpers (`read_dea()`, `read_abundances()`, `read_sites()`, …) rather than by hand.
The helpers attach the protein annotation and set the factor levels.

## Things to know

* Contrast `A_vs_B`: `diff = log2(A) − log2(B)`.
* `estimate_type = lod_imputed`: a group mean was imputed, not measured.
* Truncated total-proteome sample names (`2anc_pulseSN_1`) are repaired to `JE2anc_pulseSN_1`.
* The phospho DEAs are site-level. `site` = `<protein>_<residue><position>~<peptide>`.
* More material in the course folder, not copied here: the HTML QC/DEA reports,
  `SummarizedExperiment.rds`, `AnnData.h5ad` and `*.parquet` in each `Results_*` folder,
  plus ORA/GSEA input files.
