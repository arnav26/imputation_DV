# Imputation Pipeline

This repository contains a Snakemake pipeline for processing VCF files for imputation purposes. The pipeline involves filtering, adding tags, phasing, and imputation of reference and target datasets.

## Dependencies

1. bcftools
2. shapeit5
3. Java (for running the provided jar files)
4. Beagle (https://faculty.washington.edu/browning/beagle/beagle.22Jul22.46e.jar)
5. Conform-gt (https://faculty.washington.edu/browning/conform-gt.html)
## Configurations

The following configurations are set at the top of the Snakefile:

- `CONFORM_GT`: Path to the conform-gt jar file.
- `BEAGLE`: Path to the Beagle jar file.

## Pipeline Overview

### 1. Prepare REF

- **Filter VCF REF**: Filter the reference VCF file for variants with MAF <0.01
- **Add AN_AC_tags REF**: Add AN and AC tags to the filtered reference VCF as required by Shapeit5.
- **Shapeit5 Phasing REF**: Phase the filtered and tagged reference VCF using shapeit5 using a pedigree.
- **Convert BCF to VCF REF**: Convert the phased BCF file to VCF format for the reference.

### 2. Prepare Target

- **Conform Target**: Conform the target VCF to match the reference VCF format and allele ordering.
- **Add AN_AC_tags Target**: Add AN and AC tags to the conformed target VCF.
- **Phase Target**: Phase the conformed and tagged target VCF using shapeit5.
- **Convert BCF to VCF Target**: Convert the phased BCF file to VCF format for the target.

### 3. Imputation

- **Beagle Imputation**: Use the Beagle tool for imputing missing genotypes in the target dataset using the reference dataset.
- **Merge Imputed with Reference**: Merge the imputed target dataset with the reference dataset.

## Usage

To run the pipeline, use the following command:
