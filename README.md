# Higgs-to-Diphoton
Analysis of CMS open data regarding Higgs to Diphoton decays.


Useful, relevant tutorials: \
https://cms-opendata-workshop.github.io/workshop2021-lesson-cmssw/ \
https://opendata.cern.ch/docs/cms-getting-started-2011 
<br>
<br>

Analyzer implements several (but not all) event selection criteria from the [1], in order to limit the amount of background. 

1) In each event, two jets with the highest transverse momenta, within $|\eta| < 4.7$ are chosen. 

2) The leading and subleading of the two must have transverse momenta higher than 30 and 20 GeV respectively. 

3) Their invariant mass is required to be greater than 350 or 250 for 7 and 8 TeV datasets respectively. 

4) Furthermore, their $\eta$  separation, calculated as an absolute value of their $\eta$ difference, is required to be greater than 3.5. 



If these conditions are not met, the event is rejected. 


5) Next, two photon candidates with the highest transverse momenta, within the fiducial region $\eta$  < 2.55 excluding the barrel-end-cap  transition region, 1.44 < $\eta$ < 1.57, are chosen. 
6) Leading and subleading photons are required to have transverse momenta greater than $m_{\gamma \gamma}/3$ and $m_{\gamma \gamma}/4$ respectively. 
7) The difference between the average pseudorapidity of two jets and the pseudorapidity of a diphoton system, calculated from its total energy and momentum in z direction, must be less than 2.5. 

The invariant masses of dijets and diphotons were calculated using the equation:

$M^2 = (E_1+E_2)^2 - \left|\mathbf{p}_1 + \mathbf{p}_2\right|^2$


Additionally, diphoton invariant masses below 85 GeV were disregarded as they are much below the investigation range.
If the event passes all the required criteria, the diphoton invariant mass is added to a histogram. This event selection is passed by one good event, for approximately every 2500 to 3000 events. The histogram has 100 bins in an invariant mass range 105 to 160 GeV.

The other criteria, not implemented in this analyzer, but mentioned in [1], are:

- Difference in azimuthal angle of diphoton and dijet system must be greater than 2.6 radians.

- Photons must satisfy the ''loose'' identification criteria

From the dataset, histograms of invariant masses were obtained and saved in a ROOT file, according to the procedure described above. This data was then converted into txt format for further analysis, which was done in a Jupyter Notebook. Neighbouring histogram bins were combined, such that 25 bins remained of the initial 100, to avoid bins with low counts, as this would cause errors in uncertainty estimation and not allow weighted fit to be done properly. Background in each of the histograms was estimated through least-squares, unweighted fitting of 3rd order polynomial to the data excluding region of invariant masses between 115 GeV and 135 GeV. A Gaussian fit is then performed to the differences between data and background polynomial fit. Additionally, for each of the datasets separately, and for data from all them combined, weighted and unweighted fits of order 3 to 5 were performed. In case of weighted fit, weights were taken as inverse of vertical error. Only 2012 data produced a significant peak, with histograms from two other datasets, from 2011, being dominated by noise. 

From the 2012 run dataset with $\sqrt{s}$ = 8 TeV [2], a peak at mass of around 125 GeV was obtained for each of the 3rd - 5th order polynomial fits, with 3rd order fit giving the most significant result, and 5th, the least.


# References
[1] S. Chatrchyan et al., "Observation of a new boson at a mass of 125 GeV with the CMS experiment at the LHC," Physics Letters B, vol. 716, no. 1, pp. 30-61, Sep. 2012, doi: https://doi.org/10.1016/j.physletb.2012.08.021. <br>
[2] CMS collaboration (2022). Photon primary dataset in AOD format from Run of 2012 (/Photon/Run2012A-22Jan2013-v1/AOD). CERN Open Data Portal. DOI:10.7483/OPENDATA.CMS.2UWH.YB9E

