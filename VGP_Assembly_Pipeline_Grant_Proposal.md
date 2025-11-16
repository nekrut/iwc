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

### Core VGP Workflows (VGP0-VGP9)

| Workflow | Purpose | Key Tools | Primary Inputs | Primary Outputs | Scientific Rationale |
|----------|---------|-----------|----------------|-----------------|---------------------|
| **VGP0: Mitogenome Assembly** | Independent assembly of mitochondrial genome from HiFi reads | • MitoHiFi (reference-guided mito assembler)<br>• NCBI database queries | • Species name<br>• HiFi reads<br>• Email for NCBI<br>• Genetic code table | • Mitochondrial genome (FASTA, GenBank)<br>• Coverage plots<br>• Annotation visualization | Mitochondrial genomes serve as phylogenetic markers and can be assembled independently due to high copy number, circular structure, and conserved organization. |
| **VGP1: K-mer Profiling (Standard)** | Genome characterization through k-mer frequency analysis | • Meryl (5 operations)<br>• GenomeScope<br>• RDEval<br>• Mash | • Species name<br>• HiFi reads<br>• K-mer length (21/31)<br>• Ploidy | • Meryl k-mer database<br>• GenomeScope plots & statistics<br>• RDEval quality reports<br>• Genome size/heterozygosity estimates | K-mer analysis provides reference-free genome size estimates and reveals complexity (heterozygosity, repeats) without prior assembly. Critical for guiding assembler parameters. |
| **VGP2: K-mer Profiling (Trio)** | Enhanced genome characterization using parental data for haplotype-specific k-mer databases | • Meryl suite<br>• GenomeScope (3× profiles) | • Offspring HiFi reads<br>• Paternal Illumina reads<br>• Maternal Illumina reads<br>• K-mer length, ploidy | • 3 Meryl databases (child, paternal, maternal)<br>• 3 GenomeScope profiles with plots | Trio binning enables complete haplotype resolution via parent-specific k-mers. Valuable for highly heterozygous genomes (>1%) where standard methods struggle. |
| **VGP3: HiFi-Only Assembly** | Generate primary and alternate assemblies using only HiFi reads | • Hifiasm<br>• Cutadapt<br>• Bandage<br>• BUSCO, Compleasm, Merqury<br>• gfastats | • HiFi reads<br>• K-mer database (VGP1)<br>• GenomeScope parameters<br>• BUSCO lineage | • Primary assembly (FASTA + GFA)<br>• Alternate assembly (FASTA + GFA)<br>• QC reports (BUSCO, Merqury)<br>• Nx and size plots | HiFi reads combine long-read scaffolding with high accuracy (>99.9%), enabling diploid assembly without additional phasing data. Suitable for low-moderate heterozygosity. |
| **VGP4: HiFi + Hi-C Phasing** | Haplotype-resolved assembly using HiFi reads with Hi-C phasing | • Hifiasm (Hi-C mode)<br>• MultiQC<br>• Cutadapt<br>• BUSCO, Compleasm, Merqury<br>• gfastats | • HiFi reads<br>• Hi-C paired reads<br>• GenomeScope parameters<br>• K-mer database<br>• BUSCO lineage | • Haplotype 1 (FASTA + GFA)<br>• Haplotype 2 (FASTA + GFA)<br>• Trimmed Hi-C reads<br>• Comprehensive QC reports | Hi-C provides long-range contact information for megabase-scale phasing. Effective for moderate-high heterozygosity when parental data unavailable. |
| **VGP5: HiFi + Trio Phasing** | Haplotype-resolved assembly using parental Illumina data | • Hifiasm (trio mode)<br>• Cutadapt<br>• BUSCO, Compleasm, Merqury<br>• gfastats | • HiFi reads<br>• Paternal/maternal Illumina<br>• 3 Meryl databases (VGP2)<br>• GenomeScope parameters | • Haplotype 1 (FASTA + GFA)<br>• Haplotype 2 (FASTA + GFA)<br>• Complete QC reports | Parent-specific k-mers enable near-perfect haplotype resolution. Highest phasing accuracy, essential for highly heterozygous genomes (>1%) and haplotype-specific studies. |
| **VGP6: Purge Duplicates (Both)** | Remove artifactual duplications from both haplotypes | • purge_dups<br>• Minimap2<br>• BUSCO, Compleasm, Merqury<br>• gfastats | • Trimmed HiFi reads<br>• Primary assembly<br>• Alternate assembly<br>• K-mer database<br>• Genome size | • Purged haplotype 1 (FASTA + GFA)<br>• Purged haplotype 2 (FASTA + GFA)<br>• QC validation reports | Removes false duplications where heterozygous regions assembled separately. Transfers purged sequences between haplotypes, preserving diversity while improving quality. |
| **VGP6b: Purge Duplicates (Single)** | Purge duplications from single haplotype | • purge_dups<br>• Minimap2<br>• grep<br>• BUSCO, Compleasm, Merqury<br>• gfastats | • Assembly to purge<br>• Assembly to preserve<br>• HiFi reads<br>• K-mer database<br>• Genome size | • Purged assembly (FASTA + GFA)<br>• QC reports for both assemblies | Applied when one haplotype is high-quality and duplicates should be discarded rather than redistributed. |
| **VGP7: Bionano Scaffolding** | Scaffold contigs to chromosome-level using optical maps | • Bionano Scaffold<br>• gfastats<br>• ggplot2 | • Bionano CMAP file<br>• Genome size estimate<br>• Assembly (GFA) | • Scaffolded assembly<br>• Non-scaffolded contigs<br>• Statistics and improvement plots | Optical maps (~500kb resolution) provide orthogonal data enabling scaffolding across repetitive regions difficult to resolve by sequencing. Increases N50 and chromosome-scale contiguity. |
| **VGP8: Hi-C Scaffolding** | Generate chromosome-level scaffolds using Hi-C contact data | • YAHS<br>• BWA-MEM2<br>• Samtools suite<br>• PretextMap/Graph/Snapshot<br>• Picard, Pairtools<br>• Compleasm, MultiQC<br>• gfastats<br>(32 tools total) | • Assembly (GFA)<br>• Haplotype ID<br>• Hi-C paired reads<br>• Mapping quality threshold<br>• Restriction enzyme<br>• Genome size | • Chromosome-scale scaffolds (FASTA + GFA)<br>• Trimmed Hi-C reads<br>• Pretext contact maps (before/after)<br>• Assembly & Hi-C alignment statistics | Chromatin contacts occur more within chromosomes than between. Achieves near-complete chromosome assemblies with accurate orientation. Pretext maps validate scaffolding quality. |
| **VGP9: Decontamination** | Remove contaminants and separate mitochondrial sequences | • NCBI FCS-GX<br>• NCBI FCS-Adaptor<br>• BLAST tools<br>• Parse mito blast<br>• gfastats | • Assembly (FASTA)<br>• NCBI taxonomy ID<br>• Species binomial name<br>• Assembly/haplotype names | • Taxonomy report<br>• Contaminant sequences<br>• Mitochondrial scaffold list<br>• Adapter reports<br>• Decontaminated assembly | Contamination impacts gene prediction and evolutionary analyses. NCBI FCS tools use comprehensive databases for gold-standard detection. Essential pre-publication step. |

### Auxiliary Workflows

| Workflow | Purpose | Key Tools | Primary Inputs | Primary Outputs | Scientific Value |
|----------|---------|-----------|----------------|-----------------|------------------|
| **Plot-Nx-Size** | Visual comparison of multiple assemblies | • gfastats<br>• ggplot2<br>• datamash<br>• Text processing | Collection of FASTA assemblies (names become labels) | • Nx plot (cumulative length)<br>• Size plot (contig distributions) | Quantifies improvements from scaffolding/purging, demonstrating clear assembly progression through pipeline stages. |
| **Hi-C Contact Map for Manual Curation** | Generate Hi-C contact maps with genomic tracks for PretextView curation | • PretextMap/Graph/Snapshot<br>• Minimap2<br>• Bedtools, deeptools<br>• Seqtk<br>• gfastats<br>(17 tools total) | • Up to 2 haplotypes (FASTA)<br>• Hi-C reads<br>• Mapping quality threshold<br>• Telomere repeat<br>• PacBio reads | • Concatenated assembly<br>• Mapped Hi-C BAM<br>• Telomere/gap/coverage tracks<br>• Pretext maps (±tracks, ±filtering)<br>• Snapshot images | PretextView enables interactive breaking of mis-joins, scaffold reordering, and chromosome structure validation. Human-in-the-loop crucial for reference-quality assemblies. |

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
