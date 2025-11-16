# Vertebrate Genome Project Assembly Pipeline: A Comprehensive Computational Framework

## Executive Summary

The Vertebrate Genome Project (VGP) assembly pipeline represents a state-of-the-art computational framework for generating chromosome-level, phased genome assemblies from long-read sequencing data. This modular pipeline integrates cutting-edge bioinformatics tools within the Galaxy platform, enabling reproducible, high-quality genome assemblies that are essential for comparative genomics, evolutionary biology, and conservation research. The pipeline consists of 13 interconnected workflows that systematically process PacBio HiFi reads, Hi-C data, and optional Bionano optical maps to produce reference-quality genomes with comprehensive quality control at each stage.

## Pipeline Architecture

The VGP assembly pipeline is organized into three major functional categories:

1. **Core Assembly Workflows** (VGP0-VGP5): Initial genome characterization and contig-level assembly
2. **Assembly Refinement Workflows** (VGP6-VGP9): Purging duplicates, scaffolding, and decontamination
3. **Auxiliary Workflows**: Quality visualization and manual curation support

This modular design allows researchers to select appropriate workflows based on available sequencing data types (e.g., HiFi-only, HiFi with Hi-C, or HiFi with trio data) while maintaining consistent quality standards throughout the assembly process.

---

## Detailed Workflow Descriptions

### VGP0: Mitogenome Assembly

**Purpose**: Independent assembly of the mitochondrial genome from PacBio HiFi reads.

**Key Tools**:
- **MitoHiFi**: Specialized mitochondrial genome assembler that uses reference-guided assembly
- **NCBI database queries**: Automatic retrieval of appropriate mitochondrial reference sequences based on species taxonomy

**Inputs**: Species name, HiFi reads, email for NCBI queries, genetic code table

**Outputs**: Complete mitochondrial genome in FASTA and GenBank formats, coverage plots, and annotation visualization

**Scientific Rationale**: Mitochondrial genomes serve as important phylogenetic markers and can be assembled independently due to their high copy number, circular structure, and conserved organization across vertebrates.

---

### VGP1: K-mer Profiling for Haploid/Diploid Genomes

**Purpose**: Genome characterization through k-mer frequency analysis to estimate critical assembly parameters.

**Key Tools**:
- **Meryl** (5 operations): Fast k-mer counting and database generation with histogram generation
- **GenomeScope**: Statistical modeling of k-mer frequency distributions to infer genome properties
- **RDEval**: Read quality assessment specifically designed for PacBio HiFi data
- **Mash**: Distance-based quality control between HiFi datasets

**Analytical Outputs**:
- Genome size estimation
- Heterozygosity rate calculation
- Repeat content quantification
- Ploidy assessment
- Sequencing coverage evaluation
- Data quality metrics

**Inputs**: Species name, HiFi reads collection, k-mer length (typically 21 or 31), ploidy level

**Outputs**: Meryl k-mer database, GenomeScope plots (linear, log, transformed versions), statistical summaries, RDEval quality reports

**Scientific Rationale**: Accurate genome parameter estimation is critical for guiding assembler settings and evaluating assembly completeness. K-mer analysis provides reference-free genome size estimates and reveals genome complexity without requiring prior assembly.

---

### VGP2: K-mer Profiling for Trio Assemblies

**Purpose**: Enhanced genome characterization using parental Illumina data to partition offspring HiFi reads into haplotype-specific k-mer databases.

**Key Tools**:
- **Meryl suite**: Generates separate k-mer databases for maternal and paternal haplotypes
- **GenomeScope**: Produces three separate genome profiles (offspring and both parents)

**Inputs**: Offspring HiFi reads, paternal Illumina reads, maternal Illumina reads, k-mer length, ploidy

**Outputs**: Three Meryl databases (child, paternal haplotype, maternal haplotype), three complete GenomeScope profiles with visualization

**Scientific Rationale**: Trio binning enables complete haplotype resolution by identifying parent-specific k-mers, allowing assemblers to phase heterozygous variants with high accuracy. This approach is particularly valuable for highly heterozygous genomes where traditional assembly methods struggle.

---

### VGP3: HiFi-Only Assembly (Solo Mode)

**Purpose**: Generate primary and alternate haplotype assemblies using only PacBio HiFi long reads.

**Key Tools**:
- **Hifiasm**: State-of-the-art HiFi assembler using overlap-error-correction paradigm
- **Cutadapt**: Adapter trimming and read quality filtering
- **Bandage**: Assembly graph visualization for structural assessment
- **BUSCO**: Benchmarking Universal Single-Copy Orthologs for completeness assessment
- **Compleasm**: Alternative completeness assessment using protein domain alignment
- **Merqury**: Reference-free assembly quality evaluation using k-mers
- **gfastats**: Comprehensive assembly statistics (50 uses across all workflows)

**Quality Control Strategy**:
- Graph-based validation (Bandage visualization)
- Gene completeness assessment (BUSCO: ~3,000-13,000 genes depending on lineage)
- Protein domain completeness (Compleasm)
- K-mer accuracy and phase consistency (Merqury QV and phase blocks)
- Standard assembly metrics (N50, L50, total length, contig counts)

**Inputs**: HiFi reads, k-mer database (from VGP1), GenomeScope parameters, BUSCO lineage, assembly names

**Outputs**: Primary assembly, alternate assembly (both in FASTA and GFA formats), comprehensive QC reports, Nx and size plots

**Scientific Rationale**: HiFi reads combine long-read scaffolding capability with high base accuracy (>99.9%), enabling accurate diploid assembly without additional phasing data. This workflow is suitable for low-to-moderate heterozygosity genomes.

---

### VGP4: HiFi Assembly with Hi-C Phasing

**Purpose**: Generate haplotype-resolved assemblies using HiFi reads with Hi-C data for phasing.

**Key Tools**:
- **Hifiasm**: Integrates Hi-C contact information for haplotype phasing
- **MultiQC**: Aggregates quality metrics across multiple tools
- **Cutadapt**: Optional trimming of Hi-C reads (5bp from 5' end for Arima data)
- **BUSCO, Compleasm, Merqury**: Comprehensive quality assessment suite
- **gfastats**: Assembly statistics generation

**Hi-C Phasing Mechanism**: Hi-C reads provide long-range contact information that identifies variants occurring on the same DNA molecule (in cis), enabling accurate haplotype reconstruction over megabase scales.

**Inputs**: HiFi reads, Hi-C paired reads, GenomeScope parameters, k-mer database, BUSCO lineage, bloom filter bit setting, assembly names

**Outputs**: Haplotype 1 assembly (FASTA + GFA), Haplotype 2 assembly (FASTA + GFA), trimmed Hi-C reads, comprehensive QC reports including MultiQC aggregation

**Scientific Rationale**: Hi-C phasing is particularly effective for genomes with moderate-to-high heterozygosity where parental data is unavailable. The long-range contact information enables chromosome-scale phasing, producing more contiguous haplotype assemblies.

---

### VGP5: HiFi Assembly with Trio Phasing

**Purpose**: Generate haplotype-resolved assemblies using parental Illumina data for highly accurate phasing.

**Key Tools**:
- **Hifiasm**: Implements trio-binning algorithm using parental k-mers
- **Cutadapt**: Read preprocessing
- **BUSCO, Compleasm, Merqury**: Quality validation suite
- **gfastats**: Assembly metrics

**Trio Phasing Mechanism**: Parent-specific k-mers (hapmers) that appear in only one parent are used to partition offspring reads into maternal and paternal sets, enabling near-perfect haplotype resolution.

**Inputs**: HiFi reads, paternal/maternal Illumina reads, three Meryl databases (from VGP2), GenomeScope parameters, BUSCO lineage, haplotype names

**Outputs**: Haplotype 1 assembly (FASTA + GFA), Haplotype 2 assembly (FASTA + GFA), complete QC reports

**Scientific Rationale**: Trio-based phasing achieves the highest phasing accuracy, particularly valuable for highly heterozygous genomes (>1% heterozygosity). This approach produces the most complete haplotype separation, essential for studying structural variation and haplotype-specific biology.

---

### VGP6: Purge Duplicate Contigs (Both Haplotypes)

**Purpose**: Remove artifactual duplications resulting from unresolved haplotypic variation or overlapping contig ends.

**Key Tools**:
- **purge_dups**: Identifies and removes duplicated sequences based on read coverage and sequence similarity
- **Minimap2**: Fast all-vs-all alignment for identifying overlapping sequences
- **BUSCO, Compleasm, Merqury**: Validation that purging doesn't remove valid sequences
- **gfastats**: Updated assembly statistics

**Purging Strategy**:
1. Identify duplicates in primary assembly (haplotype 1)
2. Transfer purged sequences to alternate assembly (haplotype 2)
3. Purge duplicates from alternate assembly
4. Validate completeness is maintained

**Inputs**: Trimmed HiFi reads, primary assembly, alternate assembly, k-mer database, GenomeScope parameters, genome size estimate, BUSCO lineage

**Outputs**: Purged haplotype 1 (FASTA + GFA), purged haplotype 2 (FASTA + GFA), comprehensive QC validation

**Scientific Rationale**: Diploid assemblers can produce false duplications where heterozygous regions are assembled separately rather than being properly phased. Purging removes redundancy while preserving true haplotype diversity, improving assembly quality metrics and downstream analysis accuracy.

---

### VGP6b: Purge Duplicates (Single Haplotype)

**Purpose**: Purge duplications from a single haplotype when duplicates should be removed rather than transferred.

**Key Tools**:
- **purge_dups**: Duplicate identification and removal
- **Minimap2**: Sequence alignment
- **grep**: Text processing for filtering
- **BUSCO, Compleasm, Merqury**: Quality assessment
- **gfastats**: Statistics generation

**Use Case**: Applied when one haplotype is known to be high-quality and duplicates in the other haplotype should be discarded rather than redistributed.

**Inputs**: Assembly to purge, assembly to leave unaltered (for QC), HiFi reads, k-mer database, GenomeScope parameters, genome size, BUSCO lineage

**Outputs**: Purged assembly (FASTA + GFA), QC reports for both assemblies

---

### VGP7: Scaffolding with Bionano Optical Maps

**Purpose**: Scaffold contigs into chromosome-level assemblies using Bionano optical mapping data.

**Key Tools**:
- **Bionano Scaffold**: Hybrid scaffolder that aligns optical maps to sequence assemblies
- **gfastats**: Assembly statistics before and after scaffolding
- **ggplot2**: Visualization of scaffolding improvements

**Optical Mapping Technology**: Bionano generates genome-wide restriction maps at ~500kb resolution, bridging gaps that sequencing-based methods cannot resolve.

**Inputs**: Bionano CMAP file, genome size estimate, assembly in GFA format

**Outputs**: Scaffolded assembly, non-scaffolded contigs, statistics, Nx and size plots showing improvement

**Scientific Rationale**: Optical maps provide orthogonal information to sequencing data, enabling scaffolding across repetitive regions that are difficult to resolve with sequencing alone. This step significantly increases contiguity (N50) and chromosome-scale assembly.

---

### VGP8: Scaffolding with Hi-C Data

**Purpose**: Generate chromosome-level scaffolds using chromatin contact frequency information from Hi-C sequencing.

**Key Tools** (32 tools - most complex workflow):
- **YAHS** (Yet Another Hi-C Scaffolder): Modern Hi-C scaffolding algorithm
- **BWA-MEM2**: Fast and accurate Hi-C read alignment
- **Samtools suite**: BAM manipulation and filtering
- **PretextMap/Graph/Snapshot**: Hi-C contact map visualization
- **Picard MarkDuplicates**: PCR duplicate identification
- **Pairtools**: Hi-C pair parsing and quality filtering
- **Bamtools filter**: Advanced BAM filtering
- **Compleasm**: Completeness assessment
- **MultiQC**: Aggregate quality reporting
- **gfastats**: Assembly statistics

**Hi-C Scaffolding Principle**: Chromatin contacts occur more frequently within chromosomes than between chromosomes, creating a strong signal for ordering and orienting sequences into chromosome-scale scaffolds.

**Inputs**: Assembly (GFA), haplotype identifier, Hi-C paired reads, optional trimming, mapping quality threshold, BUSCO lineage, restriction enzyme sequence, genome size

**Outputs**: Chromosome-scale scaffolds (FASTA + GFA), trimmed Hi-C reads (if selected), assembly statistics, Hi-C alignment statistics, Compleasm report, Pretext contact maps (before/after scaffolding)

**Scientific Rationale**: Hi-C scaffolding achieves near-complete chromosome-level assemblies with accurate orientation and ordering. Pretext contact maps enable validation of scaffolding quality and identification of misassemblies, which appear as off-diagonal signals.

---

### VGP9: Assembly Decontamination

**Purpose**: Identify and remove contaminating sequences (foreign organisms, vectors, adapters) and separate mitochondrial sequences.

**Key Tools**:
- **NCBI FCS-GX** (Foreign Contamination Screening - Genome Cross-species): Detects non-target organism sequences using taxonomic classification
- **NCBI FCS-Adaptor**: Identifies sequencing adapters and artificial sequences
- **BLAST tools**: Mitochondrial sequence identification through homology search
- **Parse mito blast**: Custom parsing of mitochondrial BLAST results
- **gfastats**: Clean assembly statistics

**Contamination Sources**:
- Symbiotic or parasitic organisms
- Laboratory contaminants (bacteria, fungi)
- Vector sequences
- Adapter sequences
- Cross-species contamination

**Inputs**: Assembly (FASTA), NCBI taxonomy ID, species binomial name, assembly name, haplotype identifier, maximum mitochondrial scaffold length

**Outputs**: Taxonomy report, contaminant sequences, mitochondrial scaffold list, adapter reports (with trimming/masking recommendations), decontaminated assembly

**Scientific Rationale**: Contamination can significantly impact downstream analyses, from gene prediction to evolutionary studies. NCBI's FCS tools represent the gold-standard for contamination detection, using comprehensive databases and taxonomic awareness. This step is essential before public deposition and publication.

---

## Auxiliary Workflows

### Plot-Nx-Size: Multi-Assembly Comparison

**Purpose**: Visual comparison of multiple assemblies to track quality improvements through the pipeline.

**Key Tools**:
- **gfastats**: Assembly metrics calculation
- **ggplot2**: High-quality publication-ready plots
- **datamash**: Statistical aggregation
- **Text processing tools**: Data formatting

**Inputs**: Collection of FASTA assemblies (collection item names become plot labels)

**Outputs**: Nx plot (showing cumulative length distribution), Size plot (showing contig length distributions)

**Scientific Value**: Enables researchers to quantify improvements from scaffolding and purging steps, demonstrating clear assembly progression.

---

### Hi-C Contact Map for Manual Curation

**Purpose**: Generate detailed Hi-C contact maps with genomic feature tracks for manual curation in PretextView.

**Key Tools** (17 tools):
- **PretextMap/Graph/Snapshot**: Contact map generation and visualization
- **Minimap2**: Read alignment with optimized Hi-C parameters
- **Bedtools, deeptools**: Coverage and gap analysis
- **Seqtk**: Sequence manipulation
- **gfastats**: Assembly statistics and format conversion

**Genomic Tracks Generated**:
- PacBio read coverage (identifies low-coverage regions)
- Gap locations (shows assembly breaks)
- Telomere positions (validates chromosome termini)

**Inputs**: Up to two haplotypes (FASTA), Hi-C reads, optional trimming, mapping quality threshold, telomere repeat sequence, PacBio reads

**Outputs**: Concatenated assembly (if two haplotypes), mapped Hi-C BAM, telomere/gap/coverage tracks, Pretext maps with/without tracks (filtered and unfiltered), snapshot images

**Manual Curation Workflow**: PretextView allows researchers to interactively break mis-joins, reorder scaffolds, and validate chromosome structure based on Hi-C contact patterns. This human-in-the-loop step is crucial for achieving reference-quality assemblies.

---

## Integrated Assembly Workflows: From Reads to Chromosomes

The VGP pipeline provides three primary trajectories depending on available data types:

### Trajectory 1: HiFi-Only Assembly (Standard Pipeline)

**Workflow Sequence**:
1. **VGP1** (K-mer Profiling) → Genome characterization and parameter estimation
2. **VGP3** (HiFi-Only Assembly) → Primary and alternate contig assemblies
3. **VGP6** (Purge Duplicates) → Remove artifactual duplications
4. **VGP8** (Hi-C Scaffolding) → Chromosome-level scaffolding [requires Hi-C data]
5. **VGP9** (Decontamination) → Remove contaminants, identify mitochondrial sequences
6. **Manual Curation Workflow** → Human-guided refinement using Hi-C contact maps

**Data Requirements**: PacBio HiFi reads (primary), Hi-C reads (for scaffolding)

**Expected Outcomes**:
- Chromosome-level primary assembly
- Chromosome-level alternate assembly
- Mitochondrial genome (VGP0, run independently)
- Comprehensive quality metrics at each stage

---

### Trajectory 2: HiFi + Hi-C Phasing (Enhanced Haplotype Resolution)

**Workflow Sequence**:
1. **VGP1** (K-mer Profiling) → Genome parameter estimation
2. **VGP4** (HiFi-HiC Assembly) → Haplotype-phased assembly using Hi-C contacts
3. **VGP6** (Purge Duplicates) → Refinement of both haplotypes
4. Optional: **VGP7** (Bionano Scaffolding) → Additional scaffolding with optical maps
5. **VGP8** (Hi-C Scaffolding) → Chromosome-level scaffolding with updated Hi-C alignment
6. **VGP9** (Decontamination) → Final cleanup
7. **Manual Curation Workflow** → Interactive refinement

**Data Requirements**: PacBio HiFi reads, Hi-C reads from same individual, optional Bionano optical maps

**Expected Outcomes**:
- Two complete haplotype-resolved chromosome-level assemblies
- Superior phasing compared to Trajectory 1
- Optimal for heterozygous genomes (0.5-2% heterozygosity)

---

### Trajectory 3: HiFi + Trio (Maximum Phasing Accuracy)

**Workflow Sequence**:
1. **VGP2** (Trio K-mer Profiling) → Haplotype-specific k-mer databases
2. **VGP5** (Trio Assembly) → Parent-informed phased assembly
3. **VGP6** (Purge Duplicates) → Haplotype refinement
4. Optional: **VGP7** (Bionano Scaffolding) → Optical map scaffolding
5. **VGP8** (Hi-C Scaffolding) → Chromosome assembly [requires Hi-C]
6. **VGP9** (Decontamination) → Contamination removal
7. **Manual Curation Workflow** → Final refinement

**Data Requirements**: PacBio HiFi reads (offspring), Illumina reads (both parents), Hi-C reads (for scaffolding), optional Bionano data

**Expected Outcomes**:
- Two fully phased chromosome-level haplotypes with highest accuracy
- Near-complete separation of parental contributions
- Ideal for highly heterozygous genomes (>2% heterozygosity)
- Enables parent-of-origin studies and haplotype-specific analysis

---

## Quality Control Framework

The VGP pipeline implements multi-layered quality control at every stage:

### Assembly Completeness Assessment

**BUSCO** (Benchmarking Universal Single-Copy Orthologs):
- Searches for conserved single-copy genes expected in the lineage
- Reports: Complete (single-copy), Complete (duplicated), Fragmented, Missing
- Typical vertebrate assemblies: >95% complete BUSCOs

**Compleasm** (Complementary Assessment):
- Protein domain-based completeness using Miniprot alignment
- Often detects genes missed by BUSCO
- Provides orthogonal validation

### K-mer Based Quality

**Merqury**:
- Reference-free quality assessment using k-mer databases
- QV Score (Quality Value): base-level accuracy (target: >40, equivalent to 99.99%)
- Completeness: percentage of reliable k-mers found in assembly
- Phase blocks: continuity of haplotype phasing
- K-mer multiplicity plots: identifies collapsed/duplicated regions

### Assembly Statistics

**gfastats** (used 50+ times across workflows):
- N50/L50: Assembly contiguity metrics
- Total assembly length vs. expected genome size
- Number and size distribution of scaffolds/contigs
- GFA format validation and statistics

### Visualization-Based QC

**Bandage**: Assembly graph structure showing repeat resolution and contig connectivity

**Pretext Maps**: Hi-C contact maps revealing:
- Chromosome structure (diagonal signal)
- Misassemblies (off-diagonal patterns)
- Contamination (separate clusters)
- Inversions and translocations

**Nx and Size Plots**: Quantitative comparison showing assembly improvement through pipeline stages

---

## Technical Implementation and Reproducibility

### Galaxy Platform Integration

All workflows are implemented in Galaxy, providing:
- **Reproducibility**: Complete provenance tracking of all analyses
- **Accessibility**: Web-based interface requiring no programming skills
- **Scalability**: Can run on local servers or cloud infrastructure
- **Standardization**: Consistent tool versions and parameters
- **Interoperability**: Workflows can be exported, shared, and published

### Computational Resource Requirements

The pipeline scales from small genomes (~500 Mb) to large vertebrate genomes (3+ Gb):

**Minimum Requirements** (human-sized genome, ~3 Gb):
- CPU: 48+ cores
- RAM: 256+ GB (especially for Hifiasm and purge_dups)
- Storage: 2-3 TB for intermediate files
- Runtime: 24-72 hours depending on data volume and infrastructure

**Optimizations**:
- Bloom filters reduce memory requirements for large genomes
- Parallel processing across workflow steps
- Checkpoint-based execution allows resumption after failures

### Data Management

**Input Data Volumes** (typical vertebrate genome):
- PacBio HiFi: 30-50X coverage → 90-150 GB
- Hi-C: 50-100X coverage → 100-200 GB
- Parental Illumina (if trio): 30X each → 60-120 GB
- Bionano (if used): 100-200X coverage → 10-20 GB

**Workflow Outputs**:
- Final assemblies: 3-6 GB per haplotype
- QC reports and visualizations: 1-5 GB
- Intermediate files (can be deleted): 100-500 GB

---

## Scientific Impact and Applications

### Broad Research Applications

**Comparative Genomics**: Chromosome-level assemblies enable synteny analysis and large-scale structural variant detection across species.

**Conservation Biology**: High-quality reference genomes support population genomics studies for endangered species management.

**Evolutionary Biology**: Haplotype-resolved assemblies reveal parent-of-origin effects, allele-specific expression, and haplotype-specific structural variation.

**Medical Research**: Model organism genomes assembled with this pipeline serve as references for disease studies and functional genomics.

### Publication-Ready Outputs

The pipeline generates all analyses required for genome announcement publications:
- Assembly statistics and quality metrics
- Completeness assessments (BUSCO/Compleasm)
- Hi-C contact maps demonstrating chromosome-level scaffolding
- Comparative plots showing assembly improvement
- Contamination screening reports

### Community Standards Alignment

The VGP pipeline implements best practices defined by:
- **Earth BioGenome Project**: Standards for reference genome quality
- **NCBI GenBank**: Submission requirements including contamination screening
- **Darwin Tree of Life**: Quality thresholds for reference genomes

---

## Future Directions and Extensibility

The modular design enables straightforward integration of:
- **Ultra-long Oxford Nanopore reads**: For improved scaffolding of repetitive regions
- **PacBio HiFi-Revio data**: Higher throughput enabling deeper coverage
- **Improved scaffolding algorithms**: YAHS updates, alternative scaffolders
- **Machine learning-based QC**: Automated misassembly detection
- **Telomere-to-telomere assembly**: Gap-free chromosome assembly for model organisms

---

## Conclusion

The VGP assembly pipeline represents a comprehensive, reproducible framework for generating reference-quality vertebrate genomes. By integrating 13 specialized workflows with over 77 bioinformatics tools, it transforms raw sequencing data into chromosome-level, haplotype-resolved assemblies with extensive quality validation. The pipeline's modular design accommodates diverse data types while maintaining rigorous quality standards, making it suitable for both large-scale genome projects and individual species assemblies. Implemented in Galaxy, the pipeline ensures reproducibility and accessibility, democratizing reference genome production for the broader scientific community.

---

## References and Resources

- **VGP Tutorial**: https://training.galaxyproject.org/training-material/topics/assembly/tutorials/vgp_genome_assembly/tutorial.html
- **Galaxy Workflow Repository**: https://github.com/galaxyproject/iwc
- **Tool Citations**: Each workflow includes complete citations for all integrated tools
- **Community Support**: Galaxy community forums and VGP assembly working group

---

**Document Version**: 1.0
**Last Updated**: November 2025
**Pipeline Location**: `workflows/VGP-assembly-v2/`
