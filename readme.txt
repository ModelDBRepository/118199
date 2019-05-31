This is the readme for the TagTriC model from:

C.Clopath, L.Ziegler, E.Vasilaki, L.Buesing and W.Gerstner
Tag-Trigger-Consolidation: A Model of Early and Late
Long-Term-Potentiation and Depression
PLoS Comput Biol. 4(12): e1000248, 2008

Execute or import tagtric.py in an interactive python shell (ipython).

The function 'sim' simulates a neuron with 2 pathways of 100 input
synapses each (see U.Frey and R.Morris, Nature 385:533-536, 1997). It
takes three string arguments: 1st stimulation protocol, 2nd
stimulation, and what to plot.

- The 2 first arguments should be one of 'wtet','stet','wlfs','slfs'
  or 'nothing' which stand for weak/strong tetanus, weak/stonrg low
  frequency prpotocol and nothing, respectively.

- The 3rd one is for plotting. The string should contain: 'w1' and/or
  'w2' for the weight (decomposed in h, l and z, see the article for
  details) of the 1st and/or 2nd pathway(s); 'hl' for a two
  dimensional graph of the decomposed weights; 'p' for the protein
  production.

The function 'simoconn' simulates the experience by O'Connor
(D.O'Connor, G.Wittenberg and S.Wang, J Neurophysiol 94:1565-1573,
2005) on the frequency dependence of LTP/D. It takes a frequency (in
[Hz]) as argument.

For any question on the code contact
lorric.ziegler@epfl.ch
April 20th, 2009 - updated with corrected amplitudes for
depression/potentiation.
