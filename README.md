# VARUS: Drawing Diverse Samples from RNA-Seq Libraries
**VARUS** was originally written by Willy Bruhn as a Bachelors' thesis supervised by Mario Stanke. This repository is a copy of https://github.com/WillyBruhn/VARUS made in November 2018 and contains many bugfixes, an incremental intron database feature and an extension for using HISAT al alternative alignment program.

**VARUS** automates the selection and download of a limited number of RNA-seq reads from at NCBI's Sequence Read Archive (SRA) targeting a **sufficiently** high coverage for many genes for
the purpose of gene-finder training and genome annotation. Each iteration of the online algorithm

- selects a run to download that is expected to complement previously downloaded reads
- downloads a sample of reads ("batch") from the run with **fastq-dump**
- aligns the reads with **STAR** or **HISAT**
- evaluates the alignment


# INSTALLATION
## LINUX
Invoke the following command from the command-line in order to clone the repository: 
```sh
git clone https://github.com/MarioStanke/VARUS.git
```

**VARUS** depends on
- [samtools](http://samtools.sourceforge.net/),
- [bamtools](https://github.com/pezmaster31/bamtools), install on Ubuntu with `sudo apt-get install bamtools libbamtools-dev`
- [fastq-dump](https://ncbi.github.io/sra-tools/fastq-dump.html) and 
- [STAR](https://github.com/alexdobin/STAR) or [HISAT2](https://ccb.jhu.edu/software/hisat2) (tested with HISAT 2, version 2.0.0-beta)

Compile **VARUS** manually with
```
cd Implementation
make
``` 

### Disable NCBI Cache
By default the NCBI tool `fastq-dump` creates temporary files under `~/ncbi` of the same size as the run file from which data is downloaded, even if only a small part thereof is downloaded. Disable this caching behavior that requires probably too much hard drive space for most users with
```
mkdir -p ~/.ncbi
echo '/repository/user/cache-disabled = "true"' >> ~/.ncbi/user-settings.mkfg
```

# Getting Started

## Example
Change to directory `example` and follow the instructions in [example/README](example/README).

## Running VARUS
Copy the file `VARUSparameters.txt` from the example folder to your working directory and adjust it if necessary:

Most important parameters:

**--batchSize** specifies how many reads should be downloaded in each iteration (e.g. 50000 or 200000)

**--maxBatches** specifies how many batches should be downloaded at most

The final output is a sorted spliced alignment file (all batches together) called ***VARUS.bam***.

# References
Please cite: 
[VARUS: sampling complementary RNA reads from the sequence read archive](https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/s12859-019-3182-x). 2019; *BMC Bioinformatics*, 20:558

## VARUS at PAG2019
![Poster](docs/poster-PAG2019.png)

## Bachelor Thesis
Find the bachelor thesis of Willy Bruhn corresponding to **VARUS** in /docs/Thesis.
