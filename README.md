# Motif Scraper
Pythonic tool to search for degenerate motif matches in FASTA sequence files.

## Installation
Motif scraper is available via pip or GitHub download.
We **HIGHLY** recommend installing in a Python virtual environment.

```bash
pip install motif_scraper
```

Or user install

```bash
pip install --user motif_scraper
```

Or install from GitHub clone.

```bash
git clone https://github.com/RobersonLab/motif_scraper.git
git checkout vN.N.N # Choose highest version tag instead of vN.N.N

pip install -e .
```

## Usage

Find all sites for the CTCF motif NNDCCACYAGRKGGCASYR in GRCh38.

```bash
motif_scraper --motif NNDCCACYAGRKGGCASYR --outputFile ctcf_sites.csv --search_strand=both GRCh38.fa
```

Find CTCF sites on chromosome 1 only.

```bash
motif_scraper -r chr1 --motif NNDCCACYAGRKGGCASYR --outputFile ctcf_sites.csv --search_strand=both GRCh38.fa
```

Find CTCF sites only from position 10,000 to position 10,000,000 on chromosome 1.

```bash
motif_scraper -r chr1:10000-10000000 --motif NNDCCACYAGRKGGCASYR --outputFile ctcf_sites.csv --search_strand=both GRCh38.fa
```

Find CTCF match sites only on the top strand, using 10 processors.

```bash
motif_scraper --cores 10 --motif NNDCCACYAGRKGGCASYR --outputFile ctcf_sites.csv --search_strand=+1 GRCh38.fa
```

Search an Ensembl download of all protein coding transcript 3' UTRs for hsa-miR-10a sites on minus strand.

```bash
motif_scraper --cores 10 --motif TACCCTGTAGATCCGAATTTGTG --outputFile mir10a_sites.csv --search_strand=-1 GRCh38_3pUTRs.fa
```

Get debugging messages to troubleshoot code problems.

```bash
motif_scraper --loglevel DEBUG --cores 10 --motif TACCCTGTAGATCCGAATTTGTG --outputFile mir10a_sites.csv --search_strand=-1 GRCh38_3pUTRs.fa
```

Search for all motifs contained in a file.

```bash
motif_scraper --motif_file many_motifs.txt --outputFile many_motif_sites.csv GRCh38
```

Directly input a valid regular expression instead of a sequence motif

```bash
motif_scraper --motif N{19}CTR{3} --valid_regex GRCh38.fa
```