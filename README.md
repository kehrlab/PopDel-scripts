# PopDel-scripts
Collection of scripts used for the paper *PopDel identifies medium-size deletions jointly in tens of thousands of genomes *.

**Content**
This repository contains the scripts used for the comparisons and evaluations of the paper *PopDel calls deletions jointly in tens of thousands of genomes*. It does not contain the public data of of the [Polaris HiSeqX Diversity Cohort](https://github.com/Illumina/Polaris/wiki/HiSeqX-Diversity-Cohort), [Polaris Kids Cohort](https://github.com/Illumina/Polaris/wiki/HiSeqX-Kids-Cohort) or the [Illumina Platinum Genome NA12878](https://www.ebi.ac.uk/ena/data/view/ERR194147). They can be obtained from the respective online sources. Running the scripts might require you to adapt the paths in some of them.

 - [Deletion Data Simulation of Human Chromosome 21](#deletion-data-simulation-of-human-chromosome-21)
 - [Call Set comparison of Different Variant Callers on NA12878](#call-set-comparison-of-different-variant-callers-on-na12878)
 - [Polaris HiSeqX Diversity Cohort](#polaris-hiseqx-diversity-cohort)
 - [Polaris Kids Cohort](#polaris-kids-cohort)

 ## Deletion Data Simulation of Human Chromosome 21
 [*Simulation/deletion_data_simulation/*](https://github.com/kehrlab/PopDel-scripts/tree/master/Simulation/deletion_data_simulation)

- [simulate_deletion_haplotype.py](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/deletion_data_simulation/simulate_deletion_haplotype.py): Simulation of haplotypes.
- [simulate_deletions.py](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/deletion_data_simulation/simulate_deletions.py): Simulation of reads of the samples
- [simulate.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/deletion_data_simulation/simulate.sh): Wrapping above scripts and aligning the simulated reads to the GRCh38.
- [simulation2vcf.py](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/deletion_data_simulation/simulation2vcf.py): Writing the simulated deletions to VCF-files.

Note on random seeds: The random seed used for the simulation were 19 and 20 (for variant simulation + read-simulation of haplotype 1 and 2)  for the first sample, 20 and 21 for the second sample and so on.

 [*Simulation/delly/*](https://github.com/kehrlab/PopDel-scripts/tree/master/Simulation/delly)
 - [generateDelly.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/delly/generateDelly.sh): generates the scripts for each batch size.
 - [runDelly.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/delly/runDelly.sh): Runs the scripts generated by generateDelly.sh and measures the resource consumption.

 [*Simulation/gridss/*](https://github.com/kehrlab/PopDel-scripts/tree/master/Simulation/gridss)
 - [generateGridss.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/gridss/generateGridss.sh): Generates the scripts for each batch size.
 - [runGridss.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/gridss/runGridss.sh): Calls all the scripts generated by generateGridss.sh and measures the resource consumption.
 - [gridss.filter.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/gridss/gridss.filter.sh): Deduplicates the break ends called by gridss.

 [*Simulation/lumpy/*](https://github.com/kehrlab/PopDel-scripts/tree/master/Simulation/lumpy)
 - [generateLumpy.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/lumpy/generateLumpy.sh) and [generateSVTyper.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/lumpy/generateSVTyper.sh): Generate the necessary scripts for each batch size.
 - [runLumpy.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/lumpy/runLumpy.sh): Runs all the scripts generated by [generateLumpy.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/lumpy/generateLumpy.sh)  and measures the resource consumption.
 - [runSVTyper.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/lumpy/runSVTyper.sh): Runs all the scripts generated by [generateSVTyper.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/lumpy/generateSVTyper.sh) and measures the resource consumption.
 - [lumpy.env](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/lumpy/lumpy.env): Contains the conda environment for Lumpy and SVTyper.

 [*Simulation/podel/*](https://github.com/kehrlab/PopDel-scripts/tree/master/Simulation/popdel)
- [popdelProfile.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/popdel/popdelProfile.sh): Contains the commands for creating a profile of each bam file.
- [runPopDel.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/popdel/runPopDel.sh)): Runs all the scripts and measures the resource consumption.

[*Simulation/truth/*](https://github.com/kehrlab/PopDel-scripts/tree/master/Simulation/truth)
Contains an archive of all simulated variants for each batch size used for the evaluation. The files are the results of the deletion data simulation with above mentioned random seeds.

 [*Simulation/plots/*](https://github.com/kehrlab/PopDel-scripts/tree/master/Simulation/plots)
- [eval_bed.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/popdel/plots/eval_bed.sh): Evaluates the TP/FP/FN of all tools. Note that the range of the evaluation loop might have to be adjusted to match the batch sizes that the respective tools actually processed successfully.
- [compare_results.py](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/popdel/plots/compare_results.py): Script for alternative evaluation, based on fixed positional and size-estimate margins. can also consider genotypes of individual samples.
- [eval.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/popdel/plots/eval.sh): Wrapper for [compare_results.py](https://github.com/kehrlab/PopDel-scripts/blob/master/Simulation/popdel/plots/compare_results.py)
- simulation_plots.R: Script for generating all plots of the simulated data.

## Call Set comparison of Different Variant Callers on NA12878

[*NA12878/delly/*](https://github.com/kehrlab/PopDel-scripts/tree/master/NA12878/delly)
- [delly.environment.yml](https://github.com/kehrlab/PopDel-scripts/blob/master/NA12878/delly/delly.environment.yml): Containing the conda environment for delly.
- [Snakefile_call.delly.platinum](https://github.com/kehrlab/PopDel-scripts/blob/master/NA12878/delly/Snakefile_call.delly.platinum): Snakefile to be used with Snakemake to manage Delly's workflow.

[*NA12878/lumpy/*](https://github.com/kehrlab/PopDel-scripts/tree/master/NA12878/lumpy)
- [lumpy.environment.yaml](https://github.com/kehrlab/PopDel-scripts/blob/master/NA12878/lumpy/lumpy.environment.yml): File containing the conda environment for delly.
- [call_platinum.lumpy.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/NA12878/lumpy/call_platinum.lumpy.sh): Script for calling Lumpy and SVTyper.

[*NA12878/popdel/*](https://github.com/kehrlab/PopDel-scripts/tree/master/NA12878/popdel)
- [profile.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/NA12878/popdel/profile.sh): Contains the command used for generating the profile of NA12878.
- [platinum.profiles](https://github.com/kehrlab/PopDel-scripts/blob/master/NA12878/popdel/platinum.profiles): Contains the location of the profile generated by above command. Used as input for PopDel call.
- [Snakefile_call.popdel.platinum](https://github.com/kehrlab/PopDel-scripts/blob/master/NA12878/popdel/Snakefile_call.popdel.platinum): Snakefile to be used with Snakemake to manage the PopDel call commands for each chromosome.

[*NA12878/reference*](https://github.com/kehrlab/PopDel-scripts/tree/master/NA12878/reference)
- [contromeres.bed](https://github.com/kehrlab/PopDel-scripts/blob/master/NA12878/reference/centromeres.bed)
- [pacbio/remapped_NA12878_pacbio_deduplicated_deletions.sort.all.bed](https://github.com/kehrlab/PopDel-scripts/blob/master/NA12878/reference/pacbio/remapped_NA12878_pacbio_deduplicated_deletions.sort.all.bed):
- [personalis/remapped_NA12878_personalis_deduplicated_deletions.sort.all.bed](https://github.com/kehrlab/PopDel-scripts/blob/master/NA12878/reference/personalis/remapped_NA12878_personalis_deduplicated_deletions.sort.all.bed):

[*NA12878/plots*](https://github.com/kehrlab/PopDel-scripts/tree/master/NA12878/plots)
- [VennDiagNA12878.R](https://github.com/kehrlab/PopDel-scripts/tree/master/NA12878/plots/VennDiagNA12878.R): Creates Venn diagrams for the call sets of PopDel, Delly and Lumpy for NA12878. Works on the output of [bedtools/Snakefile.intersect](https://github.com/kehrlab/PopDel-scripts/blob/master/NA12878/plots/bedtools/Snakefile.intersect)

[*NA12878/plots/bedtools/*](https://github.com/kehrlab/PopDel-scripts/tree/master/NA12878/plots/bedtools)
- [workflow.png](https://github.com/kehrlab/PopDel-scripts/blob/master/NA12878/plots/bedtools/workflow.png): Overview of the workflow for creating the necessary overlaps for the final venn diagram.
- [workflow_generous.png](https://github.com/kehrlab/PopDel-scripts/blob/master/NA12878/plots/bedtools/workflow_generous.png): Overview of the generous workflow for creating the necessary overlaps for the venn diagram.
- [config.yaml](https://github.com/kehrlab/PopDel-scripts/blob/master/NA12878/plots/bedtools/config.yaml): Configuration of the evaluation.
- [Snakefile.intersect](https://github.com/kehrlab/PopDel-scripts/blob/master/NA12878/plots/bedtools/Snakefile.intersect): Snakefile for use with Snakemake. Manages the BED-conversion and overlap calculations via bedtools intersect.
- [Snakefile.intersect_generous](https://github.com/kehrlab/PopDel-scripts/blob/master/NA12878/plots/bedtools/Snakefile.intersect_generous): Snakefile for use with Snakemake. Manages the BED-conversion and generous overlap calculations via bedtools intersect.

## Polaris HiSeqX Diversity Cohort

[*polaris_diversity_cohort/delly/*](https://github.com/kehrlab/PopDel-scripts/tree/master/polaris_diversity_cohort/delly)
- [dag.svg](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_diversity_cohort/delly/dag.svg): Overwiev of Delly's calling workflow.
- [config.yaml](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_diversity_cohort/delly/config.yaml): Configuration for the corresponding Snakefile.
- [environment.yaml](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_diversity_cohort/delly/environment.yaml): Contains the conda environment for Delly.
- [Snakefile_polaris_delly](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_diversity_cohort/delly/Snakefile_polaris_delly): Snakefile for use with Snakemake. Manages Delly's workflow.

[*polaris_diversity_cohort/lumpy/*](https://github.com/kehrlab/PopDel-scripts/tree/master/polaris_diversity_cohort/lumpy)
- [dag.svg](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_diversity_cohort/lumpy/dag.svg): Overwiev of Lumpy's calling workflow.
- [config.yaml](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_diversity_cohort/lumpy/config.yaml): Configuration for corresponding Snakefile.
- [environment.yaml](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_diversity_cohort/lumpy/environment.yml): Contains the conda environment for Lumpy.
- [Snakefile_polaris_lumpy](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_diversity_cohort/lumpy/Snakefile_polaris_lumpy): Snakefile for use with Snakemake. Manages Lumpy's workflow.

[*polaris_diversity_cohort/popdelProfile/*](https://github.com/kehrlab/PopDel-scripts/tree/master/polaris_diversity_cohort/popdelProfile)
- [config.yaml](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_diversity_cohort/popdelProfile/config.yaml): Configuration for corresponding Snakefile.
- [Snakefile_polaris_profile](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_diversity_cohort/popdelProfile/Snakefile_polaris_profile): Snakefile for use with Snakemake. Manages the creation of the profiles for all samples.

[*polaris_diversity_cohort/popdelCall/*](https://github.com/kehrlab/PopDel-scripts/tree/master/polaris_diversity_cohort/popdelCall)
- [dag.svg](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_diversity_cohort/popdelCall/dag.svg): Overwiev of PopDel's calling workflow.
- [150.rnd.profiles](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_diversity_cohort/popdelCall/150.rnd.profiles): Shuffled list of paths to the profiles created by PopDel profile.
- [Snakefile_polaris_popdel](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_diversity_cohort/popdelCall/Snakefile_polaris_popdel): Snakefile for use with Snakemake. Applies PopDel call on each chromosome of all samples jointly.

[*polaris_diversity_cohort/plots/*](https://github.com/kehrlab/PopDel-scripts/tree/master/polaris_diversity_cohort/plots)
- [extract_calls.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_diversity_cohort/plots/extract_calls.sh): Commands for transforming the calls of the tools to the allele-count matrix required by pca_boxplot_varCount.R for the PCA.
- [extract_popdel_per_sample_variants_GT27.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_diversity_cohort/plots/extract_popdel_per_sample_variants_GT27.sh): Counts PopDel's deletions per sample.
- [ancestry.csv](https://github.com/kehrlab/PopDel-scripts/tree/master/polaris_diversity_cohort/plots/ancestry.csv): Lists the ancestry of each sample.
- [pca_boxplot_varCount.R](https://github.com/kehrlab/PopDel-scripts/tree/master/polaris_diversity_cohort/plots/pca_boxplot_varCount.R): Script for creating the PCA-plots for all tools, the box-plot for PopDel, and assessing the significance of the variant counts for PopDel.

## Polaris Kids Cohort

[*polaris_kids_cohort/delly/*](https://github.com/kehrlab/PopDel-scripts/tree/master/polaris_kids_cohort/delly)
- [dag.svg](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/delly/dag.svg): Overview of Delly's workflow.
- [enter link description here](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/delly/environment.yaml)environment.yaml: Conda environment for Delly.
- [config.yaml](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/delly/config.yaml): Contains the configuration for the corresponding Snakefile.
- [Snakefile_polaris_kids_delly](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/delly/Snakefile_polaris_kids_delly): Manages Delly's workflow.

[*polaris_kids_cohort/lumpy/*](https://github.com/kehrlab/PopDel-scripts/tree/master/polaris_kids_cohort/lumpy)
- [dag.svg](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/lumpy/dag.svg): Overview of Lumpy's workflow.
- [environment.yaml](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/lumpy/environment.yml): Conda environment for Lumpy.
- [config.yaml](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/lumpy/config.yaml): Contains the configuration for the corresponding Snakefile.
- [Snakefile_polaris_kids_lumpy](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/lumpy/Snakefile_polaris_kids_lumpy): Manages Lumpy's workflow.

[*polaris_kids_cohort/popdel/*](https://github.com/kehrlab/PopDel-scripts/tree/master/polaris_kids_cohort/popdel)
- [dag.svg](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/popdel/dag.svg): Overview of PopDel's workflow.
- [family_kids.profiles](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/popdel/family_kids.profiles): Shuffled list of paths to profiles created by PopDel profile.
- [Snakefile_polaris_kids_popdel](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/popdel/Snakefile_polaris_kids_popdel)Snakefile_polaris_kids: Snakefile for use with Snakemake. Manages PopDel's workflow.

[*polaris_kids_cohort/mendelianError/*](https://github.com/kehrlab/PopDel-scripts/tree/master/polaris_kids_cohort/mendelianError)
- [deduplicate.py](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/mendelianError/deduplicate.py): Script for removing duplicates from a VCF file.
- [kids.ped](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/mendelianError/kids.ped): Contains the pedigree information of the kids cohort. One trio <Parent1, Parent2, Child> per line.
- [mendel.py](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/mendelianError/mendel.py): Script for calculating the Mendelian error rates of the call sets.
- [mendelLoop.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/mendelianError/mendelLoop.sh): Script for applying mendel.py with varying genotype quality thresholds.
- [mendel_extended.py](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/mendelianError/mendel_extended.py): Same as mendel.py, but with additional information on error-types.
- [transmission.py](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/mendelianError/transmission.py): Script for calculating the transmission rates of the calls for varying genotype quality thresholds.

[*polaris_kids_cohort/plots/*](https://github.com/kehrlab/PopDel-scripts/tree/master/polaris_kids_cohort/plots)
 - [mendelian_plots.R](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/plots/mendelian_plots.R): Script for plotting the Mendelian error rates from the output of [mendelLoop.sh](https://github.com/kehrlab/PopDel-scripts/blob/master/polaris_kids_cohort/mendelianError/mendelLoop.sh) for all tools.
