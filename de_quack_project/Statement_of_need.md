---
title: 'de_quack: A Python library for storage, harmonization, and querying of differential gene expression results'
tags:
  - Python
  - bioinformatics
  - RNA-seq
  - differential expression
  - DuckDB
authors:
  - name: Tiger Xing
    orcid: 0009-0005-9181-7431
    affiliation: 1
affiliations:
  - name: Independent Researcher
    index: 1
date: 23 September 2026
bibliography: paper.bib
---

# Summary

de_quack is a Python library that harmonizes, stores, and queries differential gene expression (DGE) results produced by tools such as DESeq2, edgeR, and limma. It automatically reconciles inconsistent column-naming conventions across tools, resolves gene identifiers against reference nomenclature, and tracks experiment metadata for reproducibility — implemented using DuckDB and Polars for fast, file-based storage.

## Statement of Need

Differential gene expression (DGE) analysis of RNA-seq data is a widely used method for evaluating gene expression differences across conditions [@rosati2024dgereview]. 
Tools such as DESeq2 [@love2014deseq2], limma [@ritchie2015limma], and edgeR [@robinson2009edger] are well-established for this analysis, but 
their output consists of individual tables with inconsistent column naming conventions across tools, leaving data organization, cross-tool comparison, and long-term storage entirely to the user. 
This inconsistency has been identified as a real source of confusion and error, including by the developers of these tools themselves [Bioconductor support forum](https://support.bioconductor.org/p/9143655), 
and has motivated dedicated solutions within the R/Bioconductor ecosystem, such as DeeDeeExperiment [@abassi2025deedeeexperiment] and iSEEde [@ruealbrecht_isee_de], which standardize storage and retrieval of DE 
results across DESeq2, edgeR, and limma output. Broader data management frameworks such as LaminDB [@lamindb2022] address biological data storage and provenance more generally, but are not scoped specifically 
to the column-naming and identifier inconsistencies unique to DGE output. As a result, a lightweight, zero-configuration tool for harmonizing and querying DGE results specifically within the Python ecosystem 
has remained largely unaddressed.

My library, de_quack, addresses this by providing methods to efficiently ingest, query, and manipulate DGE data from various sources, such as pandas, CSV/TSV, polars, and parquet. 
Using dynamic column name mapping on ingestion, de_quack accounts for column names from Limma, edgeR, and DESeq2 along with an array of possible columns arising from custom scripts. 
On ingestion, de_quack also offers gene symbol and ensembl ID healing for the human and mouse genes through premade gene tables from the HGNC and MGI, respectively. Metadata handling 
is also stored and tracked comprehensively, with it being stored in a separate table and indexed by an experiment ID. Along with the main data storage framework, de_quack offers 
volcano plotting capabilities and a modified polars object that couples the underlying polars table with its experiment metadata by intercepting dunder methods, so the object behaves 
like a native polars DataFrame while preserving provenance through any transformation. All of these things work to make downstream DGE handling more convenient and reproducible.

de_quack is designed primarily for individual researchers, graduate students, and postdoctoral bioinformaticians managing DGE results across experiments and analysis tools, prioritizing 
zero-configuration setup and convenience over the scalability of larger, distributed data systems. Even so, the core capability of harmonizing and comparing DE results across tools and experiments 
is broadly useful wherever DGE data accumulates, regardless of project scale.

# References



