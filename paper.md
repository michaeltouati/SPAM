---
title: 'SPAM : an Open Source Code for Stopping Power of Protons and Alpha particles in Ambiant Matter'
tags:
  - Python
  - Tkinter Graphical User Interface (GUI)
  - Proton stopping power
  - Alpha particle stopping power
  - Ambient matter
  - Bragg peak
  - Continuous collisional slowing down range
authors:
  - name: Michaël J-M R TOUATI
    orcid: 0000-0001-7590-0941
    affiliation: "1"
affiliations:
 - name: Centro de Láseres Pulsados de Salamanca (CLPU), Edificio M5, Parque Cientfico, C/ Adaja 8, 37185 Villamayor, Salamanca, Spain # current affiliation
   index: 1
date: 5 August 2021
bibliography: paper.bib
---

# Summary

SPAM (Stopping Power of Protons and Alpha particles in Ambiant Matter) is a python tool using Tkinter, PIL, Numpy and Matplotlib packages that allows for visualizing, printing and saving the stopping power and/or the Bragg's peak of protons or alpha particles in ambiant matter. The code has been benchmarked against the NIST databases (https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html and https://physics.nist.gov/PhysRefData/Star/Text/ASTAR.html).

# Statement of need

`SPAM` allows for a rapid visualization, data and image generations of the stopping power and continuous slowing down range of protons and alpha particles in ambient matter using a Python Tkinter Graphical User Interface (GUI); cf \autoref{fig:spam}. The next version release will also take into account the angular scattering of protons and alpha particles in ambient matter by using a Monte-Carlo approach of a multiple binary collision model based on `@Mott:1932`, `@Moliere:1947` and/or `@Lewis:1950`.

# Mathematics

The stopping power for a proton in a material at ambient conditions 
\begin{equation}
\displaystyle \frac{ d \varepsilon }{ d s} =  {\displaystyle \left ( \frac{d \varepsilon}{ds} \right ) }_{\mathrm{ele}} + { \displaystyle \left (  \frac{d \varepsilon}{ds} \right )}_{\mathrm{nuc}}
\end{equation}
is defined as the average energy loss $d \varepsilon$ per unit path length $ds$. Due to the huge mass of atom nuclei relative to the electron mass, the proton slowing down is mainly due to Coulomb interaction of the proton with bound atomic electrons. According to Bethe theory `@Bethe:1933`, `@BetheAshkin:1953`,  the contribution of collisions with atomic electrons can be written `@ICRU49:1993`
\begin{equation}
 {\displaystyle \left ( \frac{d \varepsilon}{ds} \right ) }_{\mathrm{ele}}= 4 \pi \frac{ n_{e} e^{4}  L}{m_e v^{2}}.
\label{eq:Bethe}
\end{equation}
Here,  $e$ is the elementary charge, $m_e$ is the electron mass. $n_e$ is the atomic electron density and  $v$ is the proton velocity. The main contribution to the stopping number
\begin{equation}
L = L_0 + L_1 + L_2
\label{L}
\end{equation}
is 
\begin{equation}
L_\mathrm{0} = \displaystyle \frac{1}{2}  \ln{ \displaystyle  \left ( \displaystyle  \frac{ 2 m_{e} c^{2} {\beta^{2}} \Delta \varepsilon_{m} }{I \displaystyle \left ( 1-\beta^{2} \right ) }  \right ) } - \beta^{2}  - \frac{C}{Z} - \frac{\delta}{2} 
\label{Lo}
\end{equation}
due to the mean excitation energy $I$ of the material where the proton is propagating through. Here,
\begin{equation}
\Delta \varepsilon_\mathrm{m}=  \displaystyle  \frac{ 2 m_{e} c^{2} {\beta^{2}} }{ 1-\beta^{2} } {\displaystyle \left [ 1 + \displaystyle \frac{ 2 m_{e} }{ m_p } {\displaystyle \left (  1 - \beta^{2} \right ) }^{-  \frac{1}{2} } + {\displaystyle \left ( \displaystyle \frac{ m_{e} }{ m_p } \right )}^{2} \right ]}^{-1}
\label{W}
\end{equation}
is the largest possible energy loss by the proton in a single collision with a free electron, $m_p$ the proton mass, $\beta = v / c$ and $c$ the velocity of light in vacuum. However, as the proton kinetic energy decreases while propagating in the material, the contribution to the stopping power from interactions with bound atomic electrons in the K, L, M, ...-shells decreases and a correction term $C/Z$ must be taken into account; see `@Walske:1952` for K-shell corrections, `@Khandelwal:1968` for L-shell corrections  and `@Bichsel:1991`, `Bichsel:1992`, `Bichsel:1983` for M-shell corrections and above. Also, for relativistic proton kinetic energies, the stopping power is reduced due to the resulting electrical polarization of the medium `@Fermi:1940`, `@Sternheimer:1952`, `@Sternheimer:1982`. It is called the density effect correction because it increases with the electron density.  However, considering only non-relativistic protons, this term can be neglected in all the following.  The stopping number correction $L_1$ is the Barkas correction accounting for discrepancies between negatively and positively charged projectiles `@Barkas:1956`, `@Barkas:1963`. Finally, the second stopping number correction $L_2$ provides the valid electronic stopping power expression when the proton velocity is large compared to the velocity of bound atomic electrons `@Bloch:1933`, `@Bohr:1948`. The contribution of collisions with atomic nuclei 
\begin{equation}
{ \displaystyle \left (  \frac{d \varepsilon}{ds} \right )}_{\mathrm{nuc}} = n_{\mathrm{nuc}}  \int \Delta \varepsilon  d \sigma_{\mathrm{nuc}}
\label{Nuclear_stopping_power}
\end{equation}
to the proton slowing down is much smaller since the recoil energy 
\begin{equation}
\Delta  \varepsilon =  
\displaystyle \frac{ 4 \varepsilon \mu^2 }{ m_p m_{\mathrm{nuc}} }\sin^{2}{ \displaystyle \left ( \displaystyle \frac{ \theta}{ 2 } \right ) }  
\label{W_Nuclear}
\end{equation}
received by the target atom nucleus is small. Here,  $\mu =  m_p m_{\mathrm{nuc}} / \displaystyle \left ( m_{\mathrm{nuc}} + m_p \right )$  and $\theta$ are respectively the effective proton mass and its deflection angle in the collision center-of-mass frame, $d \Omega = 2 \pi \sin{\theta} d \theta$ and 
\begin{equation}
\label{Mott}
{\displaystyle \left ( \displaystyle \frac{ d \sigma }{ d \Omega }  \right ) }_{\mathrm{nuc}}= \displaystyle \frac{ Z^2 e^4  }{ 4 \mu^2  v^4 \sin^4{ \displaystyle \left ( \displaystyle \frac{\theta}{2} \right ) } } \displaystyle \left [ 1 - \beta^2  \sin^2{ \displaystyle \left ( \displaystyle \frac{\theta}{2} \right ) }  \right ] 
\end{equation}
is the differential cross section obtained by @Mott:1932. However, the Bethe theory (\autoref{eq:Bethe}) breaks down when the proton velocity is much lower than the orbital electron velocities. `@VarelasBiersack:1970` compiled many experimental and theoretical results `@Lindhard:1964`, `@Newton:1975`, `@AndersenZiegler:1977` and provide a fitting formula for the electronic stopping power contribution in this low velocity regime.
     
In a compound, the stopping power for a proton can be approximated by a linear combination of stopping powers in each element constituents taken separately. If we note $\rho$ the compound mass density, this so-called Bragg's additivity rule `@BraggKleeman:1905` reads 
\begin{equation}
\displaystyle \frac{d \varepsilon}{ds} = \sum_{j} \omega_{j} {\displaystyle \left ( \displaystyle \frac{d \varepsilon}{ds} \right )}_{j}.
\end{equation}
Here ${\left ( d \varepsilon / ds \right )}_j$ is the stopping power for a proton in the element constituent $j$  at ambient conditions and 
each linear coefficient $\omega_j = \mathrm{w}_j  \rho /  \rho_{j}$ depends on its fraction by weight $\mathrm{w}_j$ and its density $\rho_j$.  

# Figures

![Screenshot of SPAM GUI concerning the stopping power of non-relativistic protons in a scintillator
(vinyltoluene-based solid at room temperature) and the corresponding range of a proton with an initial kinetic
energy of 10 MeV.\label{fig:spam}](SPAM.png)

# Acknowledgements

I acknowledge the contributions from Marine Huault for having use the code in order to calibrate a HD-V2 Gafchromic films stack to proton dose response and the subsequent development of a numerical tool to characterize the 2D-resolved kinetic energy spectra of laser-generated Target Normal Sheat Acceleration protons.

# References
