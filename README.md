# editTools

> An R Package for the analysis of RNA editing data from studies utilizing whole genome sequencing and RNA-seq from the same individual(s).

## Table of Contents

1. [About](#about)
2. [Installation](#installation)
3. [Usage](#usage)
    - [Input preparation](#input-preparation)
    - [editTools run](#edittools-run)
    - [Output R object](#output-r-object)
    - [Adding RE and gene annotation](#adding-re-and-gene-annotation)
    - [Plotting](#plotting)
4. [Future plans](#future-plans)

## About

editTools was made out of the necessity to create a reproducible means to analyze RNA editing data. Candidate RNA editing detection is most often performed using scripts to analyze variant call format (VCF) data, however such scripts are rarely published leading to the possibility that variation in RNA editing results may be partly due to subtle differences in script design. This R package provides the means to analyze VCF data using standard methodology inspired from previous studies, and do so within the R framework where powerful graphical tools can be utilized to best explore the data.

From a single individual, editTools can analyze Whole genome sequencing (WGS) as well as RNAseq from any number of tissues. Input for editTools can be prepared a number of ways (allowing for differential sequencing, mapping, and variant calling techniques).

![about_diagram](./img/about_diagram.png)

An example of how I prepare input for editTools can be found [here](https://github.com/funkhou9/variant_calling_pipeline)

## Installation

> editTools contains compiled code and relies on the Rcpp package and C++11. Please note that editTools has only been tested on Mac and Linux machines and it would not be surprising if it's missing components or installation instructions needed for Windows.

If using GNU version 4.7 or later, specify that you want to install editTools using C++11 with:

```r
Sys.setenv("PKG_CXXFLAGS" = "-std=c++11")
```

Otherwise if using an older version of the GNU, use:

```r
Sys.setenv("PKG_CXXFLAGS" = "-std=c++0x")
```
Failure to use one of the above commands before attempting to install editTools will result in fatal compilation errors.

If the `devtools` package is not already installed, then install it with `install.packages("devtools")`, then install the editTools package with:

```r
devtools::install_github("funkhou9/editTools")
```

## Usage

### Input preparation

The only current requirement for editTools is that the input (resulting from the Variant Calling software of choice), must be in VCF format and contain the following header information.

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	<DNA>	<RNA1>	<RNA2>	<RNA3> ...
```

where `<DNA>` corresponds to your WGS sample, and `<RNA1>`, `<RNA2>`, `<RNA3>` ... correspond to any number of RNA samples. Currently, the INFO field must contain `DP` and `DV` tags to represent total sequencing depth and variant depth, respectively.

**VCF format is (currently, to my knowledge) unaware of strand specificity. In other words, all bases reported in a VCF file correspond to the TOP strand of the genome.** This presents a challenge when calling transcriptome variants. For RNA editing discovery, this means that one cannot distinguish an A-to-G (DNA-to-RNA) mismatch from a T-to-C mismatch, since a single cDNA read will contain both the G allele and the C allele, and knowledge of which strand that cDNA read was generated from is absent from VCF files. To workaround this and have the ability to distinguish (for example) A-to-G mismatches from T-to-C mismatches, the user can prepare **two VCF files**, one containing variants on plus-strand transcripts for each `<RNA>` sample, and one containing variants on minus-strand transcripts for each `<RNA>` sample. In each case, the same `<DNA>` sample is provided.

### editTools run

If two VCF files are prepared (see above), then identification of DNA-to-RNA mismatches can be performed with:

```r
edits <- editTools::find_edits(<plus.vcf>, <minus.vcf>, names = c("DNA", "Brain", "Liver", "Heart"))
```

- `<plus.vcf>` VCF file containing plus-strand transcripts for each RNA sample (a character string).
- `<minus.vcf>` VCF file containing minus-strand transcripts for each RNA sample (a character string).
- `names` Desired names for DNA and RNA samples (a character vector).

Without supplying additional arguments, by default editTools will scan both VCF files and identify candidate RNA edited loci according to criteria largely inspired by [Chen et al. 2014](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004274), in which a candidate RNA edited site must have:

1. A genotype that is homozygous according to 95% of the gDNA reads
2. A genotype sequencing depth of at least 10 gDNA reads
3. An "editing depth" (*the cDNA sequencing depth in support of an alternative base to what is indicated by gDNA*) of at least 5 cDNA reads

Additionally, editTools will ignore all sites containing indels in either the gDNA or any of the cDNA samples as a means to only search for single nucleotide mismatches, and will require that the gDNA sample and the cDNA sample supporting a mismatch both have a phred-scaled strand bias p-value of at least 20.

For information on all advanced options to tweak mismatch idenfication parameters, see:

```r
?editTools::find_edits
```

### Output R object

`edits` is an object of class "edit_table", and initially contains two fields: `AllSites` and `Tissues`.

* `Allsites` is a data.frame with candidate RNA editing events in the rows with the following columns:

	- `ID` A numerical ID associated with the tissue specific candidate editing event.
	- `Chr` A character representing the chromosome of the edit.
	- `Pos` A numeric representing the chromosomal position of the edit.
	- `Strand` A character representing the strand of the edited transcript ("+" or "-")
	- `Mismatch` A character representing the mismatch (eg. "AtoG")
	- `DNA_depth` A numeric representing the total DNA sample sequencing depth
	- `DNA_variant_depth` A numeric representing the DNA sample sequencing depth that is in support of a variant
	- `RNA_depth` A numeric representing the total RNA sample sequencing depth
	- `RNA_mismatch_depth` A numeric representing the RNA sample sequencing depth in support of a DNA-to-RNA mismatch
	- `RNA_edit_frac` A numeric. Simply the `RNA_mismatch_depth` / `RNA_depth`
	- `Phred_strand_bias` A numeric representing the phred-scaled probabilities of strand-bias.
	- `Ave_MQ` A numeric representing the average mapping quality at this site across all samples.
	- `Tissue` A character indicating which tissue (named with `names` argument) the event was found in.

* `Tissues` is a list with the number of elements equal to the number of tissues studied. Each element is a data.frame with mismatch types (A-to-G) in rows and the following columns:

	- `Mismatch` A character with the type of mismatch
	- `Freq` The number of events found of the specified type of mismatch
	- `Prop` The proportion of mismatches belonging to the specified type

### Adding RE and Gene Annotation

Without importing data from [RepeatMasker](http://www.repeatmasker.org/) or [Variant Effect Predictor](http://useast.ensembl.org/info/docs/tools/vep/index.html), the **edit_table** object will contain no repetitive element or gene annotation, respectively. editTools facilitates the merging of **edit_table** objects with Repeatmasker genomic datasets ([example](http://www.repeatmasker.org/species/susScr.html)) or Variant effect predictor output.

editTools implements a pseudo binary search approach to find the positions of each DNA-to-RNA mismatch within a repeatmasker dataset. In theory, this should allow for very large **edit_table** objects to be searched among large genomic repeatmasker datasets.

To add RepeatMasker annotation to the **edit_table** object `edits`, simply use:

```r
edits <- editTools::add_repeatmask(edits, <rm.out>)
```

where `<rm.out>` is RepeatMasker output from a genomic dataset, such as the one provided in the example above.

For adding gene annotation, likewise use:

```r
edits <- editTools::add_VEP(edits, <vep.out>)
```
where `<vep.out>` is VEP output generated from the positions of each mismatch.

>To generate `<vep.out>`, input for VEP can be generated with `write_vep(edits)`

>Currently VEP outputs containing the fields `Consequence`, `BIOTYPE`, `STRAND`, `SYMBOL`, `Gene`, `TREMBL`, `IMPACT`, `Amino_acids`, and `SIFT` are all required, although editTools will soon support VEP outputs containing fewer/different fields.

### Plotting

A number of plotting tools are included, such as an S3 method for `plot()` that can be used to generate a "DNA-to-RNA mismatch count histogram".

#### Some examples (after attaching editTools with `library(editTools)`):

```r
plot(edits)
```
will generate a histogram showing all DNA-to-RNA mismatch types and their counts for all tissue samples studied.

```r
plot(edits, field = "RepSites")
```
will only construct a histogram considering mismatches within repetitive elements.

Information for more plotting techniques can be found by using:

```r
edit_prop_plot(edits)
```
will produce a plot showing the distribution of editing levels among the dataset.


## Future plans

1. **Monumental changes to program design** - Incorporate new input types. The ability to process BAM files instead of VCF files could greatly enhance editTools ease of use. For instance, separating plus-strand alignments from minus-strand alignments would no longer be needed by the user. *But how would the identification of gDNA and cDNA variants be affected?*
