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
  - name: Martin Beracochea
    orcid: 0000-0003-3472-3736
    affiliation: 1
  - name: Matt Burridge
    orcid: 0000-0002-6201-2598
    affiliation: 4
  - name: Alexander Sczyrba
    orcid: 0000-0002-4405-3847
    affiliation: [3, 5]
  - name: Mahfouz Shehu
    orcid: 0009-0002-9470-0368
    affiliation: 1
  - name: Tom Tubbesing
    orcid: 0000-0003-2101-0200
    affiliation: 3
  - name: Benedikt Osterholz
    orcid: 0009-0007-0183-2799
    affiliation: [3, 5]
  - name: Anil Wipat
    orcid: 0000-0001-7310-4191
    affiliation: 4

affiliations:
  - name: EMBL-EBI
    index: 1
  - name: Institute for Translational Bioinformatics (University Hospital Tübingen) and Quantitative Biology Center (University of Tübingen)
    index: 2
  - name: Computational Metagenomics Group, Center for Biotechnology (CeBiTec), Bielefeld University
    index: 3
  - name: Newcastle University
    index: 4
  - name: IBG-5, Computational Metagenomics, Institute of Bio- and Geosciences (IBG), Research Center Juelich GmbH
    index: 5

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
Secondly, microbiome derived data typically lead to computationally expensive analyses, and so rely on the availability of high performance computing (HPC) or cloud infrastructures like the de.NBI Cloud [@citesAsAuthority:10.12688/f1000research.19013.1].
Historically, these challenges have been met within the context of individual database resources or projects [@citesForInformation:Thakur2024-ey].
These confines limit the FAIRness [@citesAsAuthority:Wilkinson2016-ig] of datasets (since they typically aren’t interlinked, directly comparable, or collectively indexed), and mean that the scope to analyse such datasets is governed by the available resources of the given project or service.
Removing these confines, by establishing a model for the federated analysis of microbiome derived data, will allow these challenges to be met by the community as a whole, and so enable major challenges like global microbial biodiversity monitoring to be met by a better connected collection of existing resources [@citesAsRecommendedReading:Waterhouse2022-wq].
More compute can be brought to bear by combining [EOSC](https://eosc.eu/) and [ELIXIR](https://elixir-europe.org/) infrastructures, [Galaxy](https://galaxyproject.org/) instances, and existing resources like [EMBL-EBI’s MGnify](https://www.ebi.ac.uk/metagenomics) [@citesAsAuthority:Richardson2023-ot], but this requires adopting a common schema for sharing analysed datasets, including their provenance. 
Such a schema can also directly contribute to the interlinking of omics layers, using research objects following a standard like RO-Crate [@citesAsAuthority:Bechhofer2013-wj; @usesMethodIn:Soiland-Reyes2022-yh] to connect linked open datasets.

We sought to address these challenges during BioHackathon Europe 2024, by designing and implementing a schema for this purpose. 
In this report, we detail our approaches, advances, and the uses of this work to allow the generation of comparable analyses on heterogeneous compute infrastructures.

Our vision of a "federated microbiome analysis service" centres on the idea that different parts of the metagenomic data flow should be able to be completed in varied ways, by varied groups, on varied infrastructure.
In practice the timeline of a metagenomic dataset may currently look like the following:

1. A sampling field expedition takes place in Location X. Metadata such as sampling locations are recorded in lab-books.
2. The sampler sends samples for DNA sequencing.
3. The sampler sends their sequences to bioinformatician collaborators for analysis with an in-house pipeline.
    1. The pipeline assembles raw whole-genome-sequencing reads into an assembled metagenome.
    2. The pipeline analyses the assembled metagenomes with a suite of taxonomic and functional analysis tools.
4. The group submits a publication detailing their work.
5. The journal requests that the data are submitted to an archive such as [ENA](https://www.ebi.ac.uk/ena) [@citesAsAuthority:Burgin2023-ds].
6. The reads data and some of the metadata are submitted to ENA.
7. An unrelated researcher wishes to know what existing metagenome data exist for Location X, and finds the raw sequencing data on ENA.
8. They request a metagenomic analysis service like MGnify to analyse the study.
9. MGnify repeat step 3, with a similar but not identical pipeline.

If a federated microbiome analysis service was sufficiently easy to opt into, then steps 4–9 could be streamlined by the original analysis data product being made publicly available (and discoverable) more directly.
Of course, existing strategies for this already exist: where publishers do not mandate that sequencing data is submitted to a dedicated repository like ENA, researchers and authors often submit their data (perhaps including the analysis products as well) to a generic digital object repository like [Zenodo](https://zenodo.org/).
This has the benefits of making the data publicly available, and resolvable from a digital object identifier ([DOI](https://www.doi.org/)).
However, it does not guarantee anything about the metadata or primary data completeness and quality, the schema, or its reusability in practice.
Likewise, the in-house pipeline used in step 3 above may not be publicly available (e.g. on [WorkflowHub](https://workflowhub.eu/)), described in any publication, or properly referenced by any analysis data products that are made available on a data repository.

Therefore in this project we have sought to create the tooling and standards necessary for groups to produce metagenomic sampling and analysis products with enough contextual metadata for them to be practically reusable.
A schematic of the envisioned federation is shown in Figure 1.
We identified several key points in the process where work was required: either in identifying and agreeing metadata terms, or in improving tooling.

![Conceptual schematic of a federated microbiome analysis service, including work done during this BioHackathon. In this scenario, metagenomics samples' and studies' metadata should be captured following a common RO-Crate profile. Pipeline processes, embedded tools, and results should be annotated so that their execution and creation provenance can be maintained. This allows, for example, for various groups to reuse assemblies whilst maintaining a chain of provenance. Analyses of metagenomics raw reads and assemblies (or metagenome-assembled genomes) should also follow a shared standard, so that downstream users, clients, and indexes can automatically understand how to index and compare these analyses from heterogeneous pipelines. Finally, metagenomic RO-Crates should make use of RO-Crate's preview features so that they can be self-rendering: for example it should be easy to view any HTML renderings of metadata and results contained within the crate.](./fig1-federated-microbiome-analysis-schematic.png)


# Methods
Any new standard or convention for publishing and ingesting metagenomic (meta)datasets must be readily adoptable by the various groups, tools, and infrastructures involved.
For this reason, we based each of our developments on existing work where possible, and considered both metadata standards (i.e. the metadata and formats that should be required and/or recommended), as well as tooling (i.e. the software to facilitate their adoption).

We therefore identified the following tracks for the BioHackathon, shown in Figure 1:

1. Define a minimum set of metadata terms for describing metagenomic datasets including studies, samples, and analyses as RO-Crates, based on existing vocabularies, terms, and minimum-information standards.
2. Develop the ability to store (as RO-Crates) the metadata associated with executions of metagenomic analysis pipeline runs across different workflow orchestration tools and languages.
3. Improve the code libraries needed to publish, ingest, and interactively explore those crates in production services.

## Rationale for using RO-Crates
RO-Crate was selected as the metadata format and packaging method for this project thanks to its basis in [JSON-LD](https://json-ld.org/) (which is widely supported with tooling in every relevant language) and its ability to reference both typed metadata entities and primary data file objects.
For example, consider an idealised structured description of a metagenomic dataset:

```yaml
# pseudocode
THIS:
    IS: "https://ebi.ac.uk/metagenomics/analysis/MGYA1.crate"
    IS_A: METAGENOMIC_CRATE
    IS_BASED_ON:
        ID: "metagenomic-sample-123"
        WHICH_IS_A: METAGENOMIC_SAMPLE
        WITH:
            location: Spain
            date: 2024
    HAS:
        METADATA_IN: 
            ID: "analysis_info.json"
            WHICH_IS_A: METAGENOMIC_ANALYSIS_DESCRIPTION
        DATA_IN:
            ID: "analysis_results.tsv"
            WHICH_IS_A: METAGENOMIC_TAXONOMIC_ANALYSIS
            WITH:
                file_path: "./data/analysis_results.tsv"
                file_type: TSV
                columns:
                    - lineage:
                        WHICH_IS_A: TAXONOMIC_ASSIGNMENT
                        reference: "https://gtdb.ecogenomic.org/"
                    - count:
                        WHICH_IS_A: OCCURRENCE_COUNT
                        source:
                            tool: "https://github.com/Ecogenomics/GTDBTk"
                            version: "2.4.0"
```

Here, terms in `UPPERCASE` signify a kind of domain-specific language in the metadata document: shared syntax that must be understood by all providers and consumers of the metagenomic crate.
Terms in `lowercase` signify dynamic values, properties or keys that may only be expected in a given narrow context: a `METAGENOMIC_SAMPLE` may be expected to come `WITH` a `location` value.

This kind of structure would allow nodes of our federated microbiome analysis service to build logic on this expected structure.
For example a service like [MGnify](https://www.ebi.ac.uk/metagenomics) [@citesAsAuthority:Richardson2023-ot] could index the existence of all available metagenomic crates:

```python
# pseudocode
def maybe_index_result(crate: Crate) -> bool:
    if crate.is_a == METAGENOMIC_CRATE:
        database.crates.add(
            id=crate.id, 
            from=crate.is_based_on.id, 
            metadata=crate.has
        )
        return True
    else:
        return False
```

A microbial biodiversity portal like [GBIF](https://www.gbif.org/) could add the taxonomic assignments from this geo-located sample to its map index:

```python
# pseudocode
def add_taxonomies_to_map(crate: Crate):
    assert crate.is_a == METAGENOMIC_CRATE
    taxonomy_file = next(
        dataset 
        for dataset in crate.has.data_in 
        if dataset.which_is_a == METAGENOMIC_TAXONOMIC_ANALYSIS
    )
    taxonomies_counts = read_csv(
        taxonomy_file.file_path,
        columns=dataset.columns.keys()
    )
    for assignment in taxonomies_counts:
    	  database.taxonomies.add(
            assignment, 
            location=crate.is_based_on.with.location
        )
```

RO-Crates are an existing mechanism for achieving this structure with only mildly more verbose syntax than this pseudocode.
In RO-Crates the metadata entities refer to (typically) Schema.org types and properties, and the file locations refer to either local paths *within* the packaged crate, or Uniform Resource Identifiers (URIs) outside the crate (i.e. on the web).
Our BioHackathon project's tracks therefore focussed on determining how to adopt existing RO-Crate standards for our use.

# Results

## Track 1: metadata standards for metagenomic crates
The basis of all metagenomic analyses is a biological sample, for example from a marine microbial environment or a human gut microbiome.
Like other bioinformatic datasets, in practice these typically belong logically to some kind of project or study.
We therefore focussed on determining the necessary metadata to describe studies and samples in RO-Crates.

### Draft minimum metadata specifications for metagenomic studies and samples
We compared the metadata standards of:

* Genomic Standards Consortium: [MIxS](https://w3id.org/mixs) (Minimum Information about any (x) Sequence (MIxS) standard) MIMS (Minimum Information about Metagenome Sequences);
* Public sequence repositories: [ENA's](https://www.ebi.ac.uk/ena) [@citesAsAuthority:Burgin2023-ds] Checklist [ERC000011](https://www.ebi.ac.uk/ena/browser/view/ERC000011) along with their Checklist implementation of MIxS standards;
* Community usage: the 100 most-frequently present metadata keys available from ENA submissions on samples that have been processed by [MGnify](https://www.ebi.ac.uk/metagenomics) [@citesAsAuthority:Richardson2023-ot].

This approach was chosen to avoid the creation of any new standard (being rooted in existing standards, implementations, and usage patterns), whilst enabling new adopters to publish compatible crates without necessarily needing to submit their data to a repository like ENA.
For example, a group collecting metagenomic datasets in a commercial context may wish to publicly publish compatible metadata associated with their sampling whilst not making publicly available the associated primary data sequences.

Our overall approach is shown in Figure 2.

![Schematic flow for how we envisage standardised metagenomic crates may be produced. Crate publishers would use various means to input, store, and curate the metadata of their studies and samples, before perhaps building an "internally useful" RO-Crate. Finally, they would convert that into an RO-Crate of the specified format by mapping to the designated Schema.org and Bioschemas terms. Where mandatory metadata are not directly available, they may use methods to infer that metadata, for example using large language models (LLMs) on less structured metadata like free-text descriptions.](./fig-x-metadata.png)

Through this process, we drafted the following metadata schema for a metagenomic study, where MUST, SHOULD, and MAY are the requirement levels for each term with meanings following [IETF RFC 2119](https://www.ietf.org/rfc/rfc2119.txt):

|        | **Field**                | **Description**                             |
| -------| ------------------------ | ------------------------------------------- |
| MUST   | `study_id`               | schema.org Identifier                       |
| MUST   | `study_title`            | Brief sequencing study description          |
| MUST   | `study_description`      | Detailed sequencing study description       |
| MUST   | `center_name`            | Submitting institution or organization name |
| MUST   | `first_created`          | Date when the study created                 |
| SHOULD | `study_alias`            | Submitter's unique identifier for the study |
| SHOULD | `broker_name`            | Organization or individual acting as broker |
| SHOULD | `status`                 | Status of the study ('public', 'private')   |
| SHOULD | `biome_identifier`       | E.g. ENVO, GOLD, ENA Tax IDs with ontology  |
| SHOULD | `samples_count`          | Number of samples in the study              |
| SHOULD | `last_update`            | Last updated date of the study              |
| SHOULD | `public_release_date`    | Intended public release date of the study   |
| SHOULD | `related_publications`   | Link to publications related to the study   |
| SHOULD | `related_geocoordinates` | Link to geocoordinate data                  |
| SHOULD | `related_studies`        | Link to other related studies               |
| MAY    | `author_name`            | Name of study contact person                |
| MAY    | `author_email`           | Email of study contact person               |

Typically, studies may have multiple alias identifiers, therefore `study_id` should be a `PropertyValue` with `sameAs` attributes for other known study accessions/IDs/URLs e.g. `ext_study_id`, `study_accession`.

We drafted the following metadata schema for a metagenomic sample:

|      | **Field**            | **Description**                                 |
| ---- | -------------------- | ----------------------------------------------- |
| MUST | `submitted_to_insdc` | submitted to INSDC                              |
| MUST | `investigation_type` | investigation type                              |
| MUST | `study_name`         | study name                                      |
| MUST | `study_id`           | id of the study to which the sample belongs     |
| MUST | `lat_lon`            | geographic location (longitude)                 |
| MUST | `geo_loc_name`       | geographic location (country and/or sea,region) |
| MUST | `collection_date`    | collection date                                 |
| MUST | `biome`              | environment (biome)                             |
| MUST | `feature`            | environment (feature)                           |
| MUST | `material`           | environment (material)                          |
| MUST | `env_package`        | environmental package                           |
| MUST | `lat_lon`            | geographic location (latitude)                  |
| MUST | `taxid`              | NCBI sample classification                      |
| MUST | `elevation`          | elevation (above ground)                        |
| MUST | `depth`              | depth (below water)                             |

### Availability of relevant Schema.org and Bioschemas terms
From these draft metadata specifications, we began searching for the availability or Schema.org and/or Bioschemas terms to adopt for each metadata field, yielding a small number of apparent matches:

* `lat_lon`, `geo_loc_name`, `related_geocoordinates` can all be of type [`https://schema.org/GeoCoordinates`](https://schema.org/GeoCoordinates), and `elevation` can be a property of it [`https://schema.org/elevation`](https://schema.org/elevation).
* `last_update` can be a property following [`https://schema.org/dateModified`](https://schema.org/dateModified).
* `material` can be a property following [`https://schema.org/material`](https://schema.org/material).
* `related_publications` can be a list of [`https://schema.org/ScholarlyArticle`s](https://schema.org/ScholarlyArticle).
* `center_name` can be an [`https://schema.org/Organization`](https://schema.org/Organization).

There appears to be a lack of appropriate Schema.org or Bioschemas terms for some of the other mandatory and recommended terms.
Notably:

* `depth` is not a possible property of [`https://schema.org/GeoCoordinates`](https://schema.org/GeoCoordinates).
* `biome` is a crucial concept in metagenomics, often following an ontology such as [GOLD](https://gold.jgi.doe.gov/ecosystem_classification) or [ENVO](https://environmentontology.com) [@citesAsAuthority:Buttigieg2016-xf], but does not appear to have a schema compatible type in existing vocabularies.

Furthermore, it would be broadly beneficial if the necessary properties to represent metagenomic studies and samples were directly available properties on Bioschemas `Studies` and `Samples` vocabulary groups.


## Track 2: metagenomic analysis pipeline execution metadata
The analysis of a metagenome typically involves executing a suite of analysis software: quality control tools, read trimming tools, taxonomic and functional analysis tools to compare sequences against reference databases, and "plumbing" scripts to link the suite of tools together [@citesAsPotentialSolution:Huson2007-ju; @citesAsPotentialSolution:Morais2022-ja; @citesAsPotentialSolution:Richardson2023-ot; @citesAsPotentialSolution:Uritskiy2018-ud].
Various workflow languages and executors are chosen to build these analysis pipelines, for example [Snakemake](https://snakemake.github.io/) [@citesAsAuthority:Molder2021-wi], Common Workflow Language (CWL) [@citesAsAuthority:Crusoe2022-sb], [Nextflow](https://www.nextflow.io) [@citesAsAuthority:Di-Tommaso2017-bx], and [Galaxy](https://galaxyproject.org/) [@citesAsAuthority:Galaxy-Community2024-ic].

To achieve interoperability between these pipelines or their results without adding further post-processing steps, metagenomic crates must be publishable from each of these workflow engines:

1. Snakemake: Snakemake is akin to a domain specific language (DSL) that adds additional syntax to the [Python language](https://www.python.org). Therefore, publishing an RO-Crate of a Snakemake pipeline execution should be possible using existing tooling, e.g. the [ro-crate-py package](https://github.com/ResearchObject/ro-crate-py) [@citesAsPotentialSolution:Chadwick2024-va].
2. CWL: Existing tooling allows a CWL-described pipeline and its execution to be published as RO-Crates [@citesAsPotentialSolution:Leo2024-wa].
3. Nextflow: Prior to this project, progress had been made towards a workflow execution RO-Crate output for Nextflow pipelines via the [`nf-prov` plugin](https://github.com/fbartusch/nf-prov/tree/workflow-run-crate).
4. Galaxy: Galaxy instances can enable support for exporting workflow execution RO-Crates [@citesAsEvidence:fair-ro-crate-in-galaxy; @citesForInformation:Hiltemann_2023].

We therefore focussed effort on Nextflow support, described next, and on potential improvements to the Python tooling which may eventually aid Snakemake support (described in Track 3 below).

### Example Nextflow pipeline

To enable quick testing of the nf-prov plugin we created a simple Nextflow pipeline based on the nf-core template [@citesAsAuthority:Ewels2020-dj; @citesAsRelated:Langer2024-pd]. The pipeline is available on GitHub at [famosab/wrrocmetatest](https://github.com/famosab/wrrocmetatest). It runs [fastp](https://github.com/OpenGene/fastp) and [megahit](https://github.com/voutcn/megahit). The README holds all necessary information to run the pipeline locally. Nf-core's tooling simplified the creation of this pipeline and so enabled us to focus our work on the nf-prov plugin.

### Process labels, ext directive and meta.yaml in nf-core modules
Process labels, ext directive and meta.yaml are three different entities which seemed fitting for out tasks of embedding tool and output descriptions of Nextflow workflows in RO-Crates. The following section explains each of those entities and shows which information can be extracted from them.

All Nextflow pipelines, but especially nf-core pipelines and pipelines created using the nf-core template, can use nf-core modules and subworkflows. These enable code reuse and modularisation of the pipeline code. Nf-core modules make use of the process labels within the module code. These labels are generalised and point towards the defined process resources limits which are defined with `conf/base.config` for example like

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
    label 'process_medium'
    // process definition here
}
```

In addition to the process label which is focussed on resource allocation, Nextflow also offers the `ext` directive which can be used to add additional information to the process. This information can be used to add metadata to the process which can be used for provenance tracking. Modules can be tagged with specific names or links to ontologies. These fields can then later be imported into the RO-Crate.
```groovy
process FASTP {
    label 'process_medium'
    ext name: 'fastp', 
        applicationCategory: '//edamontology.org/operation_0510'
    // process definition here
}
```

All nf-core modules come with a `meta.yaml` file which holds information about the tool(s) used in the module. the input and output files, a general description of the tool and the module maintainer(s). The fields relating to the used tools, the input files, and the output files are particularly relevant to the creation of the RO-Crate. The `meta.yaml` file can be used to extract the necessary information, which we added to the `nf-prov` plugin as described below.

### Embedding tool and output descriptions of Nextflow workflows in RO-Crates

A variety of metagenomics workflows exist which perform similar analytical tasks and generate comparable outputs. However, it is challenging to exchange and programmatically ingest the results from these workflows for further downstream analysis due to the lack of standardised, machine-readable descriptions. To address this gap, our approach aims to enrich metagenomics workflows with structured metadata by embedding tool and output descriptions directly into RO-Crates, specifically leveraging the [Workflow Run RO-Crate](https://w3id.org/ro/wfrun/workflow/0.5) format [@citesAsAuthority:usesMethodIn:citesAsPotentialSolution:Leo2024-wa]. This solution allows workflow outputs to be annotated with ontology terms that detail file contents, enabling interoperability and ease of reuse. It is important that existing workflows can be enhanced with minimal additions, such as specific keywords and ontology tags, without requiring any modifications to the workflow's core functionality.

After evaluating several methods, we chose to build upon a fork of the [`nf-prov` plugin](https://github.com/fbartusch/nf-prov/tree/workflow-run-crate), extending its functionality to better align with our requirements which can be found in our own [fork famosab/nf-prov](https://github.com/famosab/nf-prov/tree/workflow-run-crate) in the branch `workflow-run-crate`. This enhancement should allow for the integration of descriptions for both tools used and file contents produced within each workflow. To achieve this, we employ the Nextflow `ext` directive in each process, which specifies a unique keyword. This keyword links to metadata stored in a corresponding YAML file (e.g., `meta.yaml` for `nf-core` processes), allowing tool-specific and output-specific annotations. Workflow developers have the flexibility to include any metadata they find relevant; however, we recommend specifying at least a name, description and url for tools and tagging primary output files with ontologies that specify both file format and content. For instance, a gzipped FASTA file containing an assembly should be annotated with the format terms FASTA ([EDAM format_1929](http://edamontology.org/format_1929)) and GZIP ([EDAM format_3989](http://edamontology.org/format_3989)) as well as a data term providing information about what the file content represents ("fragment_assembly", [EDAM data_0925](https://bioportal.bioontology.org/ontologies/EDAM?p=classes&conceptid=data_0925)), adopting terms from the [EDAM ontology](https://edamontology.org) [@citesForInformation:Black2022-or].

An example YAML configuration below shows how metadata might be structured for an `assembly_process` with tools and output descriptions:

```yaml
nf_prov:
  assembly_process:
    tool:
      name: "MEGAHIT"
      description: "Single node assembler for large and complex ..."
      url: "https://github.com/voutcn/megahit"
      tool_operation: "http://edamontology.org/operation_0525"
    pigz:
      name: "pigz"
      description: "pigz is a fully functional replacement for gzip ..."
      url: "https://zlib.net/pigz/"
    outputs:
      - pattern: "k*.final.contigs.fa.gz"
        datatype: "http://edamontology.org/data_0925"
        format:
          - "http://edamontology.org/format_1929"
          - "http://edamontology.org/format_3989"
```

While substantial progress was made during the hackathon, more work remains to be done. The nf-prov plugin that our work is based on will, in some scenarios, produce invalid RO-Crates by writing nested items. We have raised this issue with the developer of the plugin. During the BioHackathon, we successfully implemented the tool description component, but output descriptions—though feasible through a similar approach—are not yet propagated to the RO-Crate. 

While our current efforts focus on Nextflow workflows, we aim to extend this approach to other workflow management systems, such as Galaxy, where a similar RO-Crate structure could provide consistent metadata across diverse platforms. Furthermore, future work would include integrating this method into existing metagenomics workflows to further promote interoperability and reuse in the metagenomics research community.

## Track 3: tooling for working with RO-Crates

### A python library for managing RO-Crates with strong type validation
One of the central challenges in handling data objects from disparate origins is ensuring that there is a shared ability to comprehend the meaning of objects and properties within them.
This is the challenge addressed by frameworks like RDF (Resource Description Framework) [@citesForInformation:Berners-Lee-gx], and by constraint languages like SHACL [@citesAsAuthority:shacl]

In the context of our metagenomic RO-Crates, this means ensuring that only entity types and properties that are themselves resolvable to a definition (e.g. a URI) are used, and ideally only terms from the definition repositories of [schema.org](https://schema.org) and [Bioschemas](https://bioschemas.org) [@citesAsAuthority:Gray2017-nx].

There have been various developments towards validation tools for RO-Crates, that ensure these standards are followed (as well as other crate mechanics, like the existence of a root dataset), for example [`rocrate-validator`](https://github.com/crs4/rocrate-validator).

However, to our knowledge all existing validation tools use *runtime* validation.
That is, the validator tool takes as input a complete crate and determines whether and how it fails to conform to the expected standards.

In other areas of software development, it is common practice to use *static* type validation: a form of static analysis where the type of variables is checked without executing the code the variables are defined in.
In practice, static type validation allows developers to see type violations as they write their code, through integrations with their Integrated Development Environments (IDEs).

In the Python ecosystem, one popular approach to introducing both static and runtime type validation is [Pydantic](https://pydantic.dev/).
Pydantic augments Python's built-in type hints (i.e., `str` for a string) with additional types (i.e., `Latitude`), and custom types (i.e., `MetagenomicSample`):

```python
from pydantic import BaseModel
from pydantic_extra_types.coordinate import Latitude, Longitude

class MetagenomicSample(BaseModel):
    lat: Latitude
    lon: Longitude
    sequence: str

my_sample = MetagenomicSample(
    lat=Latitude(52.0829),
    lon=Longitude(0.1827),
    sequence="ACTG"
)
```

Using Pydantic and a suitable IDE, a developer would know as they typed that `MetagenomicSample(lat="Europe", sequence=False)` was a type violation. 
Using additional functional validators, they could also ensure (at runtime) that `sequence` adheres to an alphabet of `[ACTG]`.

Bringing these two concepts together, we decided to adopt an existing Python library, [`pydantic2-schemaorg`](https://github.com/blurry-dev/pydantic2-schemaorg), which provides Pydantic type models for every term in Schema.org, for example using the [`schema.org/GeoCoordinates`](https://schema.org/GeoCoordinates) type to build a richer, type-checked location:

```python
from pydantic2_schemaorg.GeoCoordinates import GeoCoordinates

location = GeoCoordinates(
  longitude=14.25,
  latitude=40.808333,
  name="Sample location",
  id_="#my-sample-location"
)
```

During BioHackathon Europe 2024 we developed a proof-of-concept for a new Python package ([`pydantic-ro-crates`](https://github.com/EBI-Metagenomics/mgnify-ro-crates/tree/pydantic-ro-crates)) built on top of this.
This package introduces Pydantic models for RO-Crates and files that should be added to the RO-Crate; otherwise, the crate's metadata graph is constructed by simply appending Schema.org typed objects to the graph:

```python
roc = ROCrate()

# add the GeoCoordinates entity, defined previously, to the crate:
roc += location
```

Additional types not present in Schema.org are still supported, by building new or derived Pydantic models:

```python
class DataSetWithLocation(Dataset):
    location: GeoCoordinates

dataset = DataSetWithLocation(
    id_="./",
    name="Sample S1",
    description="My metagenomic sample",
    identifier="S1",
    location=location,
)

roc += dataset
```

The library currently supports constructing metadata graphs in this way, as well as:

1. adding files to the crate, via a Pydantic type `LocalisableFile` which stores a pointer to the file and its typed metadata;
2. rendering the RO-Crate JSON-LD;
3. packaging the RO-Crate as a `.zip`, with the metadata JSON file and other localised data files;
4. rendering an HTML preview of the metadata.

The HTML previews created by `pydantic-ro-crates` support a multi-page "website in a crate" concept.
When arbitrary HTML files are included in the crate – for example HTML reports generated by analysis pipeline tools like MultiQC – these files are automatically added to a navigation bar in the HTML preview of the crate.
This improves the usability of the human-readable crate preview, by making it easy to navigate between the rendering of the core crate metadata and renderings of other metadata like quality control information, geo-location mapping, and any data visualisations present.

In designing the new library, we imagined that some functionality (especially associated with these HTML previews) would only be desirable to certain user groups.
We therefore added a `contrib` section to the codebase, which currently contains a simple tool for generating HTML maps based on lists of `GeoCoordinates`-typed entities, using [`leaflet.js`](https://leafletjs.com/).

For example when constructing a metagenomic RO-Crate that includes metadata for a sample:

```python
from pydantic_ro_crate.contrib.mapping.render_map import render_leaflet_map
from pathlib import Path

# Use the mapping plugin to make a nice rendered map of the locations
render_leaflet_map([location], output=Path("map.html"), title=dataset.name)

# Add the map html to the crate
roc.add_localised_file(
   LocalalisableFile(
       id_="map.html",
       source_on_host=Path("map.html"),
       name="Sample map",
       description="Map of sample coordinates"
   )
)

# Package the crate as a zip: the metadata json, preview html, 
# and the included map html
roc.zip(Path("my-crate.zip"))
```

![Composite screenshot summarising the features of the newly developed `pydantic-ro-crates` library. Pythonic code is used to construct crates from strongly typed entities corresponding to Schema.org types (or those derived from them). A plugin is available to construct HTML maps. Rich multi-page HTML previews are supported.](./fig-x-pydantic-ro-crate.png)

### Crate browser
We also developed a proof of concept for a simple crate browser that can be used to view and browse the contents of an RO-Crate.
It serves the following purposes:

1. An easy and generalised way to view RO-Crate contents in a human-readable way:
   1. To provide a simple way to browse the contents of an RO-Crate using the preview.html file. This also supports a multi-page preview. i.e if the crate contains multiple HTML files, the browser will display a list of all the HTML files in the crate and allow the user to navigate between them.
2. A simple way to convert any RO-Crate to a self contained browsable website
   1. By leveraging on the new `pydantic_ro_crate` library, we can easily convert any RO-Crate to a self-contained website that can be easily shared with others. By wrapping this functionality of the library in a simple web server, we can provide a simple way to convert any RO-Crate to a self-contained website that can be easily shared with others.
   2. This can be deployed as a hosted service, that can be used to convert any RO-Crate to a self-contained website that can be easily shared with others.
3. The RO-Crate Browser can also be used as an independent JavaScript library that can be embedded in any web page to provide a simple way to view the contents of an RO-Crate. We've published 3 variants of this on NPM
   1. React Library: [react-ro-crate-browser](https://www.npmjs.com/package/react-ro-crate-browser)
   2. Vue Library: [vue-ro-crate-browser](https://www.npmjs.com/package/vue-ro-crate-browser)
   3. Web Component: [ro-crate-browser-web-component](https://www.npmjs.com/package/ro-crate-browser-component)
   
Figure 4 shows the architecture of the RO-Crate Browser:

![The architecture of the `RO-Crate Browser` ](./fig-x-ro-crate-browser-architecture.png)


### Converting an RO-Crate to a Browsable self-contained website:
The Crate browser built on previous work done by David Lopez on the [ro-crate-zip-explorer repository](https://github.com/davelopez/ro-crate-zip-explorer/tree/main/examples/vue/ro-crate-zip-vue]). 
This was a tool initially designed to list the contents of an RO-Crate zip file.
As mentioned above, the RO-Crate browser builds on this to give a more interactive and user-friendly way to view the contents of an RO-Crate zip file as a self-contained website.
The two diagrams below demonstrate how it works:

![An image describing the RO-Crate conversion feature of the `RO-Crate Browser` ](./fig-x-converting-ro-crates.png)

### Browsing an RO-Crate as a self-contained website:
![An image describing the browsing feature of the `RO-Crate Browser` ](./fig-x-ro-crate-browsing.png)


# Discussion
Enabling interoperability between multiple metagenomic analysis pipelines and their outputs can reduce the duplicated effort and resource utilisation needed for microbiome research.
In this BioHackathon project, we conceptualised this as a "federated microbiome analysis service" in which multiple research groups and data services could share intermediate- and end- data products from their pipelines: metadata descriptions, metagenomic assemblies, metagenome-assembled genomes, and taxonomic/functional analysis profiles.
We identified and worked on metadata standards and analysis tooling toward this vision, adopting RO-Crates [@citesAsAuthority:Bechhofer2013-wj; @usesMethodIn:Soiland-Reyes2022-yh] as the standard to support it.

## RO-Crate types, linking, and hierarchies
RO-Crates are a generic approach to packaging research objects and their metadata, whether those objects be primary data files like metagenomic sequencing reads, downstream data files like a metagenome's taxonomic profile, or a workflow codebase like an analysis pipeline.
Interlinking crates of these various types is done by referencing each by its URI and type: a metagenomic analysis crate could have a root [`schema.org/Dataset`](https://schema.org/Dataset) that [`isBasedOn`](https://schema.org/isBasedOn) a metagenomic sample crate at `doi.org/10.1000/example/sample`.
[Workflow Run RO-Crates](https://w3id.org/ro/wfrun/workflow/0.5) [@citesAsAuthority:usesMethodIn:citesAsPotentialSolution:Leo2024-wa] uses a [`schema.org/CreateAction`](https://schema.org/CreateAction) to semantically link a workflow (which is itself a `Dataset`) to its output `Dataset`.
Formalised workflow datasets (e.g. the released codebase of a pipeline) can be made available through [WorkflowHub](https://workflowhub.eu/), a public repository for workflows.
For example [https://doi.org/10.48546/workflowhub.workflow.384.3](https://doi.org/10.48546/workflowhub.workflow.384.3) points to the workflow "metaGOflow: A workflow for marine Genomic Observatories' data analysis", and a [Workflow RO-Crate](https://w3id.org/workflowhub/workflow-ro-crate/) can be rendered for this workflow.
For primary datasets, some repositories support RO-Crate export of data/metadata objects, notably [Dataverse](https://dataverse.org/) through external plugins [@citesAsPotentialSolution:Bloemen2024-jb].

In the metagenomic context, primary datasets typically follow a hierarchy similar to study – samples – assemblies – analyses.
Whilst objects at each layer of this hierarchy could be rendered as individual RO-Crates following the RO-Crate schema and adopting recommended terms for metagenomics, in practice a single RO-Crate for a study, that was itself hierarchically structured, may facilitate easier interoperability.
Future work on recommended metagenomic crate schemas would therefore benefit from support for multiple "entry-points" in RO-Crates, which is currently under discussion for version 2 of the RO-Crate specification (e.g. through Project 19 of BioHackathon Europe 2024).

## Adoption of metadata standards
Beyond the tooling needed to export RO-Crates, two distinct areas of metadata annotation are needed for adopters of a proposed metagenomic crate:

1. primary metadata annotation, for example geospatial, sampling time, and sample processing metadata;
2. computational tool annotation, for example EDAM terms for the software tools, versions, purposes, and output files of each meaningful step in a pipeline.

In practice these are likely to be the biggest blockers to adoption of this format, due to the initial and marginal work required to annotate existing and new datasets and workflows.
This can be eased by taking a pragmatic approach to the metadata recommendations: introducing no mandatory metadata fields beyond those either needed by common tools (like sequencing technology type to select an appropriate metagenome assembler tool), required for repository submission, or required by the community's adopted minimum standards.
Designing standards to be achievable with low effort will also help, similar to how nf-core's requirement that modules and submodules must emit a channel of version numbers/labels for each tool used requires only adding code like `fastqc --version` to the pipeline.

## Repository support in Europe
In Europe, and especially within the context of [ELIXIR](https://elixir-europe.org/) resources, arguably the three most fundamental repositories for the proposed metagenomic crates are:

1. [ENA](https://www.ebi.ac.uk/ena) / [BioSamples](https://www.ebi.ac.uk/biosamples) [@citesAsAuthority:Burgin2023-ds], the canonical deposition repository for primary metagenome sequences and their sampling/study metadata;
2. [EuropePMC](https://europepmc.org/) [@citesAsAuthority:Rosonovski2024-sg], a central index of life-science publications and pre-prints;
3. [WorkflowHub](https://workflowhub.eu/), a repository of computational workflows that is the de-facto choice of the European bioinformatics community.

Support for schema.org/Bioschemas terms in ENA/Biosamples would be advantageous both for data deposition, and for data retrieval.
Ideally, every BioSample and ENA Project/Study would have an RO-Crate compatible endpoint describing the entity's submitted and curated metadata, and new sequencing datasets could be submitted by uploading an RO-Crate containing the sequencing data and metadata descriptors.
This is currently not possible, although support for submissions via ISA-JSON [@citesAsAuthority:Sansone2012-ud] is being developed by the [Multi-omics Adapter for Repository Submissions (MARS)](https://github.com/elixir-europe/MARS/) project, and the similarity in structure and vocabulary between ISA-JSON and schema.org terms mean this is a promising avenue for future work.
Europe PMC does not currently providate an RO-Crate endpoint for publications, however the construction of a [`schema.org/ScholarlyArticle`](https://schema.org/ScholarlyArticle) object for a publication with a Europe PMC URL is trivial.
As previously discussed, WorkflowHub already supports rendering a deposited workflow as an RO-Crate.
In conclusion, these fundamental repositories either have existing support for RO-Crates or support a reasonable approach to building an integration.

## Future work
This BioHackathon Europe 2024 project leads naturally to several areas of further work:

1. feature completeness of the `pydantic-ro-crates` and various RO-Crate browser code libraries, beyond proofs-of-concept;
2. publishing of a draft metagenomic crate specification for studies, samples, and analyses;
3. publishing specifications for missing schema types needed to represent minimum metagenomic metadata standards;
4. tooling to convert metagenomic crates to ISA-JSON for repository submission;
5. RO-Crate export support on existing metagenomic databases like MGnify;
6. publishing the updated `nf-prov` plugin for tracing nextflow workflow.

Within this project we have investigated and validated the technical feasibility of a federated microbiome analysis services running across myriad research datasets, workflows, and infrastructures.
We have also identified areas where future effort would push this federated service toward reality.


## Acknowledgements
We thank ELIXIR, the research infrastructure for Life-science data, and the organisers of BioHackathon Europe 2024 for delivering the BioHackathon and funding the travel costs of ABR and AS.

We are grateful to Eli Chadwick and Nick Duty for fruitful discussions about the future directions of RO-Crate and Bioschemas, respectively.

## References
