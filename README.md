# Higgs-to-Diphoton
Analysis of CMS open data regarding Higgs to Diphoton decays.


Useful, relevant tutorials: \
https://cms-opendata-workshop.github.io/workshop2021-lesson-cmssw/ \
https://opendata.cern.ch/docs/cms-getting-started-2011


Analyzer implements several (but not all) event selection criteria from the \cite{16}, in order to limit the amount of background. 
\begin{enumerate}
\item In each event, two jets with the highest transverse momenta, within $|\eta| < 4.7$ are chosen.

\item The leading and subleading of the two must have transverse momenta higher than 30 and 20 GeV respectively. 

\item Their invariant mass is required to be greater than 350 or 250 for 7 and 8 TeV datasets respectively. 

\item Furthermore, their $\eta$  separation, calculated as an absolute value of their $\eta$ difference, is required to be greater than 3.5. 

\end{enumerate}
