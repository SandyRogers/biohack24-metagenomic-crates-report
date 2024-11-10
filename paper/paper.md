---
title: 'Enhancing multi-omic analyses through a federated microbiome analysis service'
title_short: 'BioHackEU24 #21: federated microbiome analyses'
tags:
  - metagenomics
  - rocrates
  - federated
authors:
  - name: Alexander B. Rogers
    orcid: 0000-0002-4283-6135
    affiliation: 1
  - name: Famke Bäuerle
    orcid: 0000-0003-1387-0251
    affiliation: 2
  - name: Tom Tubbesing
    orcid: 0000-0003-2101-0200
    affiliation: 3
  - name: Add Yourselves
    orcid: 0000-0000-0000-0000
    affiliation: 4
affiliations:
  - name: EMBL-EBI
    index: 1
  - name: University of Tuebingen
    index: 2
  - name: Bielefeld University
    index: 3
  - name: Fourth Affiliation
    index: 4
date: 8 November 2024
cito-bibliography: paper.bib
event: BH24EU
biohackathon_name: "BioHackathon Europe 2024"
biohackathon_url:   "https://biohackathon-europe.org/"
biohackathon_location: "Barcelona, Spain, 2024"
group: Project 21
# URL to project git repo --- should contain the actual paper.md:
git_url: https://github.com/SandyRogers/biohack24-metagenomic-crates-report
# This is the short authors description that is used at the
# bottom of the generated paper (typically the first two authors):
authors_short: Rogers \emph{et al.}
---

# Introduction
Multi-omics datasets are an increasingly prevalent and necessary resource for achieving scientific advances in microbial ecosystem research. 
However, they present twin challenges to research infrastructures: firstly the utility of multi-omics datasets relies entirely on interoperability of omics layers, i.e. on formalised data linking. 
Secondly, microbiome derived data typically lead to computationally expensive analyses, i.e. on the availability of powerful compute infrastructures. 
Historically, these challenges have been met within the context of individual database resources or projects. 
These confines limit the FAIRness of datasets (since they typically aren’t interlinked, directly comparable, or collectively indexed), and mean the scope to analyse such datasets is governed by the available resources of the given project or service. 
Removing these confines, by establishing a model for the federated analysis of microbiome derived data, will allow these challenges to be met by the community as a whole. 
More compute can be brought to bear by combining EOSC and ELIXIR infrastructures, Galaxy instances, and existing resources like EMBL-EBI’s MGnify, but this requires adopting a common schema for sharing analysed datasets, including their provenance. 
Such a schema can also directly contribute to the interlinking of omics layers, using research objects to connect linked open datasets.

We sought to address these challenges during BioHackathon Europe 2024, by designing and implementing a schema for this purpose. 
In this report, we detail our approaches, advances, and the uses of this work to allow the generation of comparable analyses on heterogeneous compute infrastructures.

Our vision of a ``federated microbiome analysis service'' centres on the idea that different parts of the metagenomic data flow should be able to be completed in varied ways, by varied groups, on varied infrastructure.
In practice, a metagenomic project's data flow may follow a timeline like the following:

1. A sampling field expedition takes place in Location X. Metadata such as sampling locations are recorded in lab-books.
2. The sampler sends samples for DNA sequencing.
3. The sampler sends their sequences to bioinformatician collaborators for analysis with an in-house pipeline
  a. The pipeline assembled raw whole-genome-sequencing reads into an assembled metagenome
  b. The pipeline analyses the assembled metagenomes with a suite of taxonomic and functional analysis tools
4. The group submits a publication detailing their work
5. The journal requests that the data are submitted to an archive such as ENA
6. The reads data and some of the metadata are submitted to ENA
7. An unrelated researcher wishes to know what existing metagenome data exist for Location X, and finds the raw sequencing data on ENA
8. They request a metagenomic analysis service like MGnify to analyse the study
9. MGnify repeat step 3, with a similar but not identical pipeline

If a federated microbiome analysis service was sufficiently easy to opt into, then steps 4–9 could be streamlined by the original analysis data product being made publicly available more directly.
Of course, existing strategies for this already exist: where publishers do not mandate that sequencing data is submitted to a dedicated repository like ENA, researchers and authors often submit their data (perhaps including the analysis products as well) to a generic digital object repository like Zenodo.
This has the benefits of making the data publicly available, and resolvable from a DOI.
However, it does not guarantee anything about the metadata or primary data quality, schema, or reusability in practice.

Therefore in this project we have sought to create the tooling and standards necessary for groups to produce metagenomic sampling and analysis products with enough contextual metadata for them to be practically reusable.
A schematic of the envisioned federation is shown in Figure 1.
We identified several key points in the process where work was required: either in identifying and agreeing metadata terms, or on tooling.

![Conceptual schematic of a federated microbiome analysis service, including work done during this BioHackathon. In this scenario, metagenomics samples' (and studies') metadata should be captured following a common RO-Crate profile. Pipeline processes, embedded tools and results should be annotated so that their execution and creation provenance can be maintained. This allows, for example, for various groups to reuse assemblies whilst maintaining a chain of provenance. Analyses of metagenomics raw reads and assemblies (or metagenome-assembled genomes) should also follow a shared standard, so that downstream users, clients, and indexes can automatically understand how to index and compare these analyses from heterogeneous pipelines. Finally, metagenomic RO-Crates should make use of RO-Crate's preview features so that they can be self-rendering - for example it should be easy to view any HTML renderings of metadata and results contained within the crate.](./fig1-federated-microbiome-analysis-schematic.png)



# Exemplary Nextflow pipeline

To enable quick testing of the nf-prov plugin we created a simple Nextflow pipeline based on the nf-core template [cite nf-core paper here]. The pipeline is available at [this GitHub repository](https://github.com/famosab/wrrocmetatest). It runs [fastp](https://github.com/OpenGene/fastp) and [megahit](https://github.com/voutcn/megahit). The README holds all necessary information to run the pipeline locally.

### Process labels in nf-core pipelines
All nf-core pipelines and pipelines created using the nf-core template make use of the process labels within the module code. These labels are generalized and point towards the defined process resources limits which are defined with `conf/base.config` for example like
```groovy
withLabel:process_medium {
    cpus   = { 6     * task.attempt }
    memory = { 36.GB * task.attempt }
    time   = { 8.h   * task.attempt }
}
```
which can be used in the module code like
```groovy
process FASTP {
    tag "$meta.id"
    label 'process_medium'
    // process definition here
}
```

<!-- 
# Formatting

This document use Markdown and you can look at [this tutorial](https://www.markdowntutorial.com/).

## Subsection level 2

Please keep sections to a maximum of only two levels.

## Tables and figures

Tables can be added in the following way, though alternatives are possible:

Table: Note that table caption is automatically numbered and should be
given before the table itself.

| Header 1 | Header 2 |
| -------- | -------- |
| item 1 | item 2 |
| item 3 | item 4 |

A figure is added with:

![Caption for BioHackrXiv logo figure](./biohackrxiv.png)

# Other main section on your manuscript level 1

Lists can be added with:

1. Item 1
2. Item 2

# Citation Typing Ontology annotation

You can use [CiTO](http://purl.org/spar/cito/2018-02-12) annotations, as explained in [this BioHackathon Europe 2021 write up](https://raw.githubusercontent.com/biohackrxiv/bhxiv-metadata/main/doc/elixir_biohackathon2021/paper.md) and [this CiTO Pilot](https://www.biomedcentral.com/collections/cito).
Using this template, you can cite an article and indicate _why_ you cite that article, for instance DisGeNET-RDF [@citesAsAuthority:Queralt2016].

The syntax in Markdown is as follows: a single intention annotation looks like
`[@usesMethodIn:Krewinkel2017]`; two or more intentions are separated
with colons, like `[@extends:discusses:Nielsen2017Scholia]`. When you cite two
different articles, you use this syntax: `[@citesAsDataSource:Ammar2022ETL; @citesAsDataSource:Arend2022BioHackEU22]`.

Possible CiTO typing annotation include:

* citesAsDataSource: when you point the reader to a source of data which may explain a claim
* usesDataFrom: when you reuse somehow (and elaborate on) the data in the cited entity
* usesMethodIn
* citesAsAuthority
* citesAsEvidence
* citesAsPotentialSolution
* citesAsRecommendedReading
* citesAsRelated
* citesAsSourceDocument
* citesForInformation
* confirms
* documents
* providesDataFor
* obtainsSupportFrom
* discusses
* extends
* agreesWith
* disagreesWith
* updates
* citation: generic citation -->


# Results

## Embedding Tool and Output Descriptions of Nextflow Workflows in Research Object Crates

A variety of metagenomics workflows exist which perform similar analytical tasks and generate comparable outputs. However, it is challenging to exchange and programmatically ingest the results from these workflows for further downstream analysis due to the lack of standardized, machine-readable descriptions. To address this gap, our approach aims to enrich metagenomics workflows with structured metadata by embedding tool and output descriptions directly into Research Object (RO) Crates, specifically leveraging the Workflow Run RO Crate format (cite: https://doi.org/10.1371/journal.pone.0309210). This solution allows workflow outputs to be annotated with ontology terms that detail file contents, enabling interoperability and ease of reuse. It is important that existing workflows can be enhanced with minimal additions, such as specific keywords and ontology tags, without requiring any modifications to the workflow's core functionality.

After evaluating several methods, we chose to build upon a fork of the `nf-prov` plugin (https://github.com/fbartusch/nf-prov/tree/workflow-run-crate), extending its functionality to better align with our requirements. This enhancement should allow for the integration of descriptions for both tools used and file contents produced within each workflow. To achieve this, we employ the Nextflow `ext` directive in each process, which specifies a unique keyword. This keyword links to metadata stored in a corresponding YAML file (e.g., `meta.yaml` for `nf-core` processes), allowing tool-specific and output-specific annotations. Workflow developers have the flexibility to include any metadata they find relevant; however, we recommend specifying at least a name, description and url for tools and tagging primary output files with ontologies that specify both file format and content. For instance, a gzipped FASTA file containing an assembly should be annotated with the format terms FASTA (EDAM format_1929) and GZIP (EDAM format_3989) as well as a data term providing information about what the file content represents ("fragment_assembly", EDAM data_0925)

An example YAML configuration below shows how metadata might be structured for an `assembly_process` with tools and output descriptions:

```
nf_prov:
  assembly_process:
    tool:
      name: "MEGAHIT"
      description: "Single node assembler for large and complex metagenomics NGS reads ..."
      url: "https://github.com/voutcn/megahit"
      tool_operation: "http://edamontology.org/operation_0525"
    pigz:
      name: "pigz"
      description: "pigz is a fully functional replacement for gzip that exploits ..."
      url: "https://zlib.net/pigz/"
    outputs:
      - pattern: "k*.final.contigs.fa.gz"
        datatype: "http://edamontology.org/data_0925"
        format: ["http://edamontology.org/format_1929", "http://edamontology.org/format_3989"]
```

While substantial progress was made during the hackathon, more work remains to be done. The nf-prov plugin that our work is based on will, in some scenarios, produce invalid RO crates by writing nested items. We have raised this issue with the developer of the plugin. During the biohackathon, we successfully implemented the tool description component, but output descriptions—though feasible through a similar approach—are not yet propagated to the RO crate. 

While our current efforts focus on Nextflow workflows, we aim to extend this approach to other workflow management systems, such as Galaxy, where a similar RO Crate structure could provide consistent metadata across diverse platforms. Furthermore, future work would include integrating this method into existing metagenomics workflows to further promote interoperability and reuse in the metagenomics research community.

# Discussion

...

## Acknowledgements

...

## References
