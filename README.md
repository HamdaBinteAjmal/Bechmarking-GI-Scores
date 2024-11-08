Folder "InputData" contains raw counts and associated files for each of the following CRISPR CDKO studies:<br/>
Dede,M. et al. (2020) Multiplex enCas12a screens detect functional buffering among paralogs otherwise masked in monogenic Cas9 knockout screens. Genome Biol., 21.<br/>
Parrish,P.C.R. et al. (2021) Discovery of synthetic lethal and tumor suppressor paralog pairs in the human genome. Cell Rep., 36, 109597.<br/>
Gonatopoulos-Pournatzis,T. et al. (2020) Genetic interaction mapping and exon-resolution functional genomics with a hybrid Cas9--Cas12a platform. Nat. Biotechnol., 38, 638–648.<br/>
Ito,T. et al. (2021) Paralog knockout profiling identifies DUSP4 and DUSP6 as a digenic dependence in MAPK pathway-driven cancers. Nat. Genet., 53, 1664–1672.<br/>
Thompson,N.A. et al. (2021) Combinatorial CRISPR screen identifies fitness effects of gene paralogues. Nat. Commun., 12.<br/>

The remaining folders are code scripts to apply each scoring method to each dataset: Parrish Score, zdLFC, Gemini, Orthrus and Horlbeck score. <br/>
For Gemini and Orthrus, their respective R packages have been used to calculate GI scores. <br/>
zdLFC and Parrish score has been implemented by adapting the code provided by the Dede et al. and Parrish et al. respectively.<br/>
For Horlbeck score, we have used the API provided by the SLKB project (https://slkb.osubmi.org/) <br/>
