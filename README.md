# Fontan_2025
Code to reproduce Fontan 2025

To reproduce most analyses you will need a directory with the sorted and indexed bam files.
because of the compexity of working with multiple reference

we recommand that you put your bam file into separate directory based on the reference used
Moreover we followed this convention:
  trimmed_genotype_sample_lane_VSreference.bam
  with genotype (sim / sim_mau/ mau)
  
If you face any difficulty please open an issue on this github repo.


We used Rust, python3.6+, and R 4+ in this study.

### You will need to replace path by the actual path on your computer.

## 1 - Coverage Plot
We used an in house tool with the last version available at https://github.com/rLannes/coverage_rust_bam.
We provides the Code of the version used for this paper in this directory.
you need to install it first, please refer to the readme file inside the directory "coverage_rust_bam".
Once you have it installed, run the coverage notebook "CoverPlotSI.ipynb". you will need the aligned bam.


## 2 - OmniSplice Plot
you will need to run omnisplice (https://github.com/rLannes/OmniSplice; https://www.biorxiv.org/content/10.1101/2025.04.06.647416v1) with default parameter and  "--read-to-write soft-clipped  --spliced-def 9
--unspliced-def 10 11 12 13" option
You will need to run omnisplice with D.mauritiana Y genes annotated AND the refeseq simulans (for the X and autosomes).

#### 1 - To reproduce the plot of the different splicing defect:
run the jupyter notebook:
  "OmniSplice_table_plotSI.ipynb"

#### 2 - To make the plot that show splicig defect with intron length run:
run SplicingEfficiencyPlotSI.ipynb

#### 3 - backsplicing
run the backsplicing module of omnisplice on the Currated Y gene annotation and the Dryad genome. 
To do so you will need bowtie2 and the bowtie2 reference of the genome you ran omnisplice.
In this study we ran it on the Y linked genes in the Ching-Ho assembly.


## 3 - DE


## 4 - Dot Plot
To make the dot plot we implemented a custom programs.
you'll need rust installed with maturin ( pip install maturin pip install numpy)

Please follow those preparatory step:

1- Compile aligner:
cd GlobalAlignerDNA
cargo build --release

2- compile and build python wheel.
cd ../pairAlnWindow_dotplot
maturin build --release 

3- add maturin wheel to your python (optionaly you cna use virtual env)
python -m pip install -U ./target/wheels/*.wheel

4- finally run the dot plot notebook
dotPlotPlot.ipynb

## 5 - exon coverage 
To reproduce the exon coverage plot.


### Align the read and sort and index the bam file:

first you need to align the read as described in the manuscript, to both D. mauritiana and D. simulans dryad genome.
sort an index the bam files


### gtf formatting

Then you need to merge the gtf file with our manually currated annotation.

for D. Mauritiana:
with "dmau_gene.gtf" the dryad gtf
and "Currated_Mau_Ygenes_formatted.gtf" the manually annotated genes
run:
bedtools intersect -v  -a dmau_gene.gtf  -b  Currated_Mau_Ygenes_formatted.gtf > NolargeY.gtf
cat dmau_gene.gtf NolargeY.gtf > dmau.gtf
Do the same for D. Simulans

#### extract intron and exons
in the notebook : intron_exon_coverageSI.ipynb
ran the part get intron exon gtf

#### gtf sorting



then sort the bed file as described in bedtools documentation

then we need to reorder the gtf to match the chromosomes order of the bam files so we can use the sorted option of bedtools


1- generate the genome file from the bam file header: ( do this for both the D. mau and D. sim bam files)
	samtools view -H <bam_file> | grep "@SQ" | sed -r 's/^@SQ\tSN://' | sed 's/LN://' > genome_file_dmau_dryad.txt


For this you need to run the function  match_bed_to_genome_file()
 # in the notebook part match bed to genome

This function take three argument in this order:
1- the path to the sorted gtf file
2- the path to the matching genome file
3- the path to the output gtf file


example:
match_bed_to_genome_file( "dmau_intron.sorted.bed",
                         "genome_file_dmau_dryad.txt",
                          "dmau_intron.bamsorted.bed")

##### bedtools coverage
Finally you are ready to run bedtools coverage
you need to ran it for both D.mau and D.sim

example for D.mau
bedtools coverage -sorted -g genome_file_dmau_dryad.txt -a dmau_exon.bamsorted.bed -b dmau_dryad_SRR22548176.sorted.bam > coverage_dmauExon_dmau_SRR22548176

the coverage file "coverage_dmauExon_dmau_SRR22548176" contains all the informatio we need

you can now run the plot part of the notebook.




