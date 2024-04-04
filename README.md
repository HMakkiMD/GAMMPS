# GAMMPS
GAMMPS contains several Python modules by which one can calculate the electronic structure properties of semiconducting polymers from their chemical drawings.

**Step 1** calculates force field parameters for a given repeat unit made from a sequence of conjugated monomers and sidechains all are specified in _parameters.py_ file. 
It consists of five different modules:
* _OligomerBuilder.py_ 
* _ChargeTorsionInput.py_ 
* _Oeff.py_ 
* _TorsionCorrection.py_ 
* _FF_test.py_.
 
The input files for each conjugated monomer in the repeat unit have to be provided:
* Optimised coordinate files in (i) _xyz_ and (ii) _mol_ formats 
* an _atp_ file containing the LJ parameters of all atoms (with the same sequence of lines as in the coordinate files).

**Step 2** constructs the coordinate (_xyz_ format) and force filed files for the repeat unit structure with/out sidechains and polymers with any specified length in _parameters.py_. It contains three modules: 
* _RU_builder.py_
* _SC_RU_builder.py_
* _PO_builder.py_.

**Step 3** creates input files for electronic structure calculations for "melt" and "soup" simulations as explained in (DOI). It contains three modules: 
* _DOSindex.py_
*  _DOSinput_melt.py_
* _DOSinput_soup.py_

**Step 4** performs DFT calculations (by _Gaussian 16_) to obtain the Density of State (DOS) and Localisation Length (LL) of the polymers. It contains one module:
* _QC_calculation.py_

The input files generated in **Step 3** should be provided in a folder (i.e., _input_files_) and placed in the same folder as the code. Input variables should also be defined in _input_variables.inp_ file. This module parallelises multiple DFT calculations to be efficiently performed on a single computer node with 40 CPU cores in the current version.
