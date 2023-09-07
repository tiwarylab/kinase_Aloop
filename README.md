# kinase_Aloop

## Purpose:
Please note that this is intended to be a code and data repository for a specific set of systems (i.e. protein kinases)
For a view of the full pipeline excluding multiple walker learning and static bias, please use the repository for [AF2RAVE](https://github.com/tiwarylab/alphafold2rave/)

## Motivation
[Arxiv link will be added soon.]
This protocol is essentially a methodology that combines two schools of thought: structure prediction, and enhanced sampling to preserve thermodynamics.
In particular, we aim to make the process of sampling kinase activation loop conformational diversity smoother and faster, and demonstrate this idea on DDR1 wild type and 3 mutants.

## Components

Movies, bias files, and trajectories can be found [in this Drive](https://drive.google.com/drive/folders/1XyCVlFqbUwnPY3J1f1P36W4E4e-AeOiI)

In folder mdscripts
* `openmm_utils.py` - functions required to run this instance of our methodology
* `basicmd.py` - to run basic unbiased MD
* `equilprotocol.py` - to run equilibration protocol as in the paper
* `distancemetad.py` - to run a metadynamics simulation using OP computed with SPIB with distances as input CVs
* `kinaseCVs.py` - to extract the CVs we use

To run the code, biases need to be downloaded from the drive linked above and run with the code in distancemetad.py

## Citation

Please cite the following reference if using this protocol with or without the provided code:

<<add arxiv link soon>>

* "AlphaFold2-RAVE: From sequence to Boltzmann ensemble"
Bodhi P. Vani, Akashnathan Aranganathan, Dedi Wang, Pratyush Tiwary
J. Chem. Theory Comput. 2023; doi: https://doi.org/10.1021/acs.jctc.3c00290

* "State predictive information bottleneck", Dedi Wang and Pratyush Tiwary, J. Chem. Phys. 154, 134111 (2021) https://doi.org/10.1063/5.0038198



## License

MIT License

Copyright (c) 2023 Tiwary Research Group

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
