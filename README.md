# Electrode models

In ERT, electrodes (typically consisting of well-conducting steel) are mostly considered point electrodes and modelled as such as their real extent is small compared to the electode distance. However, sometimes this is not the case and some sort of discretization including the electrode conductivity will better describe the reality. Generally there are three different ways:
- conductive cell model (CCM): the electrode volume is discretized by cells of good conductivity 
- complete electrode model (CEM): the electrode surface is discretized and coupled through contact impedances
- shunt electrode model (SEM): the electrode is described by a number nodes and connecting edges implying a shunt resistance

![Electrode models after Ronczka et al. (2015a)](SEM_CEM_CCM_sketch.svg)

We published some papers on this subject, first introducing the CEM (Rücker & Günther, 2011), followed by a study of SEM in comparison with CEM and CCM (Ronczka et al., 2015a). Finally, field-scale examples using SEM electrodes for shallow (Ronczka et al., 2015b) or deep (Günther et al., 2015) saltwater intrusion were presented. Later on, the CEM was used to model and invert data from SAMOS, a subsurface electrode chain using ring electrodes (Ronczka et al., 2020, 2025). 

The underlying codes and data were later ported to pyGIMLi and collected here.

## References
- Rücker, C. & Günther, T. (2011): The simulation of Finite ERT electrodes using the complete electrode model. *Geophysics* 76(4), F227-238, [doi:10.1190/1.3581356](https://doi.org/10.1190/1.3581356).
- Ronczka, M., Rücker, C. & Günther, T. (2015a): Numerical study of long electrode electric resistivity tomography – Accuracy, sensitivity and resolution. *Geophysics* 80(6), E317-328, [doi:10.1190/geo2014-0551.1](https://doi.org/10.1190/geo2014-0551.1).
- Ronczka, M., Voss, T. & Günther, T. (2015b): Cost-efficient imaging and monitoring of saltwater in a shallow aquifer by using long-electrode ERT. *J. of Appl. Geoph.* 122, 202-209, [doi:10.1016/j.jappgeo.2015.08.014](https:/doi.org/10.1016/j.jappgeo.2015.08.014).
- Günther, T., Ronczka, M. & Voß, T. (2015): Saltwater Monitoring Using Long-Electrode ERT. in: *Liebscher, A. & Münch, U. (Eds.): Geological Storage of CO2 – Long Term Security Aspects*, Springer Int. Publ., 167-182, [doi:10.1007/978-3-319-13930-2_8](https://doi.org/10.1007/978-3-319-13930-2_8)
- Ronczka, M., Günther, T., Grinat, M. & Wiederhold, H. (2020): Monitoring freshwater-saltwater
interfaces with SAMOS - installation effects on data and inversion. *Near Surface Geophysics* 18(4),
369-383, [doi:10.1002/nsg.12115](https://doi.org/10.1002/nsg.12115)
- Ronczka, M., Günther, T., Grinat, M., Siemon, B., Scheihing, K., & Müller-Petke, M. (2025).
Multi-scale electrical resistivity imaging and long-term monitoring as a tool for groundwater
management. *Near Surface Geophysics*, [doi:10.1002/nsg.70006](https://doi.org/10.1002/nsg.70006)