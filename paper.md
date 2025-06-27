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
  - name: Cristina Rea
    orcid: 0000-0002-9948-2649
    affiliation: 1
  - name: MIT PSFC Disruptions Group
    affiliation: 1
affiliations:
 - name: MIT Plasma Science and Fusion Center, Cambridge MA, USA
   index: 1
date: 27 June 2025
bibliography: paper.bib
---

# Summary

Magnetically-confined fusion experiments routinely operate under a wide variety of engineering parameters in order to gain invaluable insight into fusion plasmas with the purpose of understanding and then harnessing their intrinsic power for energy production.
Such exploration of a wide parameter space sometimes results in unexpected and rapid loss of confinement of the plasma discharge, events which are generically known as 'disruptions'.
Disruptions represent a significant danger to both modern experimental machines and, above all, future reactor-relevant devices.
Therefore preventing disruptions, detecting them, and avoiding them are features of paramount importance for any plasma control system.
Given the sheer number of available diagnostic systems, and possible plasma modeling tools, artificial intelligence and machine learning (AI/ML) models are ideal candidates for heavy-duty numerical computation.
Fast and agile numerical frameworks for database preparation and preprocessing are necessary for letting researchers focus on novel algorithms and benchmark different architectures and models.

# Statement of need

`DisruptionPy` [@trevisan_2025; @rea_2024; @trevisan_2024; @wei_2024] is an open-source physics-based scientific framework for disruption analysis of fusion plasmas, designed with the explicit purpose of streamlining database preparation of experimental fusion data to allow efficient AI/ML workflows.
`DisruptionPy` originated as an institutional effort from the Plasma Science and Fusion Center within the Massachusetts Institute of Technology (MIT PSFC) to create a shared and validated set of feature-extraction routines, and evolved into an open-source scientific framework in order to aid disruption scientists everywhere.
`DisruptionPy` natively supports efficiently extracting data from `MDSplus` [@stillerman_1997], the leading open-source storage back-end for most fusion experiments, and enables scientists to carry out complicated Python-based computations at scale across entire experimental databases.
`DisruptionPy` relies on established numerical libraries, e.g. `NumPy`, `SciPy`, `Pandas`, `Xarray`, to allow effortless manipulation of either raw or pre-processed data into complicated feature-extraction workflows for database generation.

The heterogeneous set of scripts from which `DisruptionPy` was developed led to several high-profile scientific publications [@hu_2021; @keith_2024; @maris_2024; @montes_2019; @rea_2018; @rea_2018; @rea_2019; @rea_2020; @spangher_2025; @tinguely_2019; @zhu_2020; @zhu_2021; @zhu_2023].
`DisruptionPy` itself is now the basis for the scientific work of the entire Disruptions Group at MIT PSFC and will undoubtedly lead to further high-impact results in the near future.

# Acknowledgements

The most recent revamp of DisruptionPy was partially supported by DOE FES under Award DE-SC0024368, "Open and FAIR Fusion for Machine Learning Applications".

# References

