# Data for 'Selective Control of Surface Spin Current in topological materials based on Pyrite-type OsX2 (X = Se, Te) Crystals'

This repository provides the minimal data sets and codes associated with the publication [**Selective Control of Surface Spin Current in topological materials based on Pyrite-type OsX2 (X = Se, Te) Crystals**](https://www.nature.com/articles/s41535-019-0186-8), consisted of following:

* Relaxed coordinates file (VASP POSCAR format) for OsSe2 (bulk and slab)
* Band structure raw output (band_non_soc.dat and band_soc.dat)
* wannier90_hr.dat (Wannier TB model) file generated from wannier90 (projected on to d orbitals of Os and p orbitals of X). You need [wannier90](http://www.wannier.org/) and/or [WannierTools](https://wanniertools.com) to do process this file to generate results such as topological invariants and surface staes.
* A customized Fortran90 code ported to WannierTools (2.3.1) to generate planar surface states specturm as shown in Figure 5(c).

More information can be obtained upon reasonable requests from the corresponding authors of the publication.

## Licence and Use

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is partly licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

These data are released under licence [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/legalcode) except for the customized Fortran 90 code. The Fortran 90 code follows the licence rule as described in [GNU General Public License V3](https://www.gnu.org/licenses/gpl-3.0.en.html). If you use this code, please kindly cite the WannierTools publication in [Here](https://doi.org/10.1016/j.cpc.2017.09.033) and [**Selective Control of Surface Spin Current in topological materials based on Pyrite-type OsX2 (X = Se, Te) Crystals**](https://www.nature.com/articles/s41535-019-0186-8).


