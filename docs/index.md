# DisruptionPy

A Python package for plasma disruption analysis and prediction. 

## Background
A key element of plasma control systems (PCS) in tokamak reactors is the prediction and avoidance of disruptions, sudden losses of the thermal and magnetic energy stored within the plasma that can occur when tokamaks operate near regions of plasma instability or because of system malfunctions. The energy released during  disruptions can cause severe damage to plasma-facing components, limiting experimental operation or even the device lifetime. This poses a serious challenge to next-step fusion experiments such as SPARC, which will have to operate near some of the limits of plasma stability to achieve its intended performance and will do so at for long and frequent intervals. Previous work has shown the promise of machine-learning (ML) algorithms for disruption prediction in both DIII-D and EAST -- the Experimental Advanced Superconducting Tokamak in China -- PCS. This is also due to the fact that fusion science currently lacks first-principle, theoretical solutions to fully predict and avoid disruptions. 

DisruptionPy is an open-source python package for training, updating, and evaluating algorithms for disruption prediction and avoidance that can be applied to DIII-D and EAST PCSs and inform SPARC disruption management strategies. It is being a developed and maintained by MIT EECS Meng student as part of his thesis "A machine-learning driven framework for plasma disruption detection in tokamaks."

## Project layout
    disruption_py # Source code
    docs # Mkdocs generated documentation
    iris_requirements # requirements.txt for D3D iris cluster
    matlab # Original matlab scripts
    ml_log # Tracking thesis experiments 
    notebooks # Example notebooks for analysis and visualization
    requirements # Requirements for locally installed environment
    scripts # Scripts for various disruption_py supported workflows
## References and Resources 
- [Design Document](https://probable-argument-b7b.notion.site/Workflow-Design-Document-a04529032bda4a999f42e75182a43258)
- [Thesis Proposal](https://www.overleaf.com/read/xyhqcgvzssqb)
- Roadmap(TODO) 