---
title: 'DisruptionPy: An open-source physics-based scientific framework for disruption analysis of fusion plasmas'
tags:
  - python
  - plasma physics
  - disruptions
  - nuclear fusion
  - magnetic confinement
authors:
  - name: Gregorio L. Trevisan
    orcid: 0000-0003-4920-6109
    affiliation: 1
    corresponding: true
  - name: Yumou Wei
    orcid: 0000-0002-8868-0017
    affiliation: 1
  - name: Amos M. Decker
    orcid: 0009-0004-3728-3259
    affiliation: 1
  - name: Joshua Lorincz
    orcid: 0009-0003-0288-9501
    affiliation: 1
  - name: Samuel L. Jackson
    orcid: 0000-0001-5301-5095
    affiliation: 2
  - name: Cristina Rea
    orcid: 0000-0002-9948-2649
    affiliation: 1
  - name: Robert S. Granetz
    orcid: 0000-0002-6560-1881
    affiliation: 1
  - name: MIT PSFC Disruptions Group
    affiliation: 1
affiliations:
 - name: MIT Plasma Science and Fusion Center, Cambridge MA, USA
   index: 1
 - name: U.K. Atomic Energy Authority, Culham Centre for Fusion Energy, Culham Science Centre, Abingdon, U.K.
   index: 2
date: 27 June 2025
bibliography: paper.bib
---

# Summary

`DisruptionPy` [@trevisan_2026; @trevisan_2024; @rea_2024; @wei_2024] is an open-source physics-based scientific framework for disruption analysis of fusion plasmas, designed with the explicit purpose of streamlining database preparation of experimental fusion data to allow efficient AI/ML workflows.

`DisruptionPy` originated as an institutional effort from the Plasma Science and Fusion Center within the Massachusetts Institute of Technology (MIT PSFC) to create a shared and validated set of feature-extraction routines, and evolved into an open-source scientific framework in order to aid disruption scientists everywhere.
`DisruptionPy` natively supports efficiently extracting data from `MDSplus` [@stillerman_1997; @stillerman_2025], the leading open-source storage back-end for most fusion experiments, and enables scientists to carry out complicated Python-based computations at scale across entire experimental databases.
`DisruptionPy` also supports extracting data using Xarray, which is interoperable with the open MAST [@jackson_2024; @jackson_2025] dataset, therefore enabling researchers to easily access and analyze historical MAST data without the need to participate in a collaboration agreement.
`DisruptionPy` relies on established numerical libraries, e.g. `NumPy`, `SciPy`, `Pandas`, `Xarray`, to allow effortless manipulation of either raw or pre-processed data into complicated feature-extraction workflows for database generation.

The heterogeneous set of scripts from which `DisruptionPy` was developed led to several high-profile scientific publications [@hu_2021; @keith_2024; @maris_2024; @montes_2019; @rea_2018; @rea_2018; @rea_2019; @rea_2020; @spangher_2025; @tinguely_2019; @zhu_2020; @zhu_2021; @zhu_2023].
`DisruptionPy` itself is now the basis for the scientific work of the entire Disruptions Group at MIT PSFC and will undoubtedly lead to further high-impact results in the near future.

# Statement of need

Magnetically-confined fusion experiments routinely operate under a wide variety of engineering parameters in order to gain invaluable insight into fusion plasmas with the purpose of understanding and then harnessing their intrinsic power for energy production.
Such exploration of a wide parameter space sometimes results in unexpected and rapid loss of confinement of the plasma discharge, events which are generically known as 'disruptions'.
Disruptions represent a significant danger to both modern experimental machines and, above all, future reactor-relevant devices.
Therefore preventing disruptions, detecting them, and avoiding them are features of paramount importance for any plasma control system.
Given the sheer number of available diagnostic systems, and possible plasma modeling tools, artificial intelligence and machine learning (AI/ML) models are ideal candidates for heavy-duty numerical computation.
Fast and agile numerical frameworks for database preparation and preprocessing are necessary for letting researchers focus on novel algorithms and benchmark different architectures and models.

As the Fusion Community prepares for the upcoming burning plasma devices, the multiple existing data repositories already face numerous interoperability challenges.
Previous community reporting [@humphreys_2020] identified the need to improve several aspects of existing platforms, ranging from hardware and technology to software, including development of optimized ML-ready workflows for scientific discovery.
Therein, the authors highlighted the current different data access systems, the various data storage formats, and a lack of adequately-labeled data as main challenges that need to be addressed by the research community.

In such context, the open-source development of `DisruptionPy` satisfies the crucial need for shared and validated data-processing workflows, and the framework's helpfulness will only grow as more experimental devices relax their requirements for data access and evolve towards open data and FAIR (Findable, Accessible, Interoperable, Reusable) principles [@wilkinson_2016].

Additional example of similar frameworks for experimental data retrieval and database preparation are `TokSearch` [@sammuli_2018] and `DEFUSE` [@pau_2023].
The `TokSearch` library [@sammuli_2018] was developed to efficiently query, process, and analyze experimental data from DIII-D for ML applications.
It leverages a distributed file format to increase throughput and a dedicated API to transfer data from MDSplus and export it in Parquet format.
`TokSearch` appears to be established only for DIII-D workflows.
`DEFUSE` [@pau_2023], the Disruption and Event analysis framework for FUSion Experiments developed, implements an interface layer to access the data from different fusion experiments through MDSplus and HDF5.
Source data, diagnostics, machine descriptions, and data-processing schemes are defined in interoperable data libraries in JSON format within a data abstraction layer.
DEFUSE has been applied to several devices, however the framework has not been open-sourced yet.

# Acknowledgements

The original MATLAB development efforts were funded by the MIT Plasma Science and Fusion Center Magnetic Confinement Fusion Experiment Research and Related Activities DE-SC0014264.
The most recent revamp of DisruptionPy was partially supported by DOE FES under Award DE-SC0024368, "Open and FAIR Fusion for Machine Learning Applications".

# References

