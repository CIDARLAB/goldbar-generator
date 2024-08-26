# GOLDBAR-Generator
This software tool generates the GOLDBAR (Grammars for cOmbinatoriaL Design and Bio-Assembly Revision) and Categories for Constellation based on known principles of Genetic Circuit Design.

#
## Rules
P - Promoter Road Blocking: Road-blocking promoters must be alone or in front of regular promoters

L - Leaky Terminators: Leaky Terminators can only be at end of circuit

O - Orthogonality: Certain parts can not be in the same circuit

I - Part-junction Interference: Certain parts can not be next to each other 

R - Do not Reuse Certain Parts: Certain parts will not be reused in the same sequence

#
## Cloning the REPO
```
git clone https://github.com/CIDARLAB/goldbar-generator.git
cd goldbar-generator/
```

#
## Installing Python Packages
Install using pip
```
pip install pandas
pip install numpy
```

#
## Using the GOLDBAR Generator
Run
```
goldbar-generator.py
```
To Run in Terminal
```
python goldbar-generator.py
```

Follow the prompts for which Part Library you would like to use then the number of transcriptional units and rules.

#
### Library/Sample Inputs
Can be found in the [library](/library/) folder. This includes the part library, not orthogonal library, and part junction interference library

#
### Example Output
Can be found in the [output](/output/) folder. This example is based on the libraries in [library](/library/) folder and uses P L O I R rules.

#
### Try Constellation
To see the Design space graph and enumerate designs use [Constellation](https://github.com/CIDARLAB/constellation-js.git).

#
### Try Knox
To save and perform further graph analysis, use [Knox](https://github.com/CIDARLAB/knox.git). You can output your design from Constellation as SBOL and import this SBOL file into Knox.

#
### References
Brophy, J., Voigt, C. Principles of genetic circuit design. Nat Methods 11, 508–520 (2014). https://doi.org/10.1038/nmeth.2926

# 
### Manuscript

Nicholas Roehner, James Roberts, Andrei Lapets, Dany Gould, Vidya Akavoor, Lucy Qin, D. Benjamin Gordon, Christopher Voigt, and Douglas Densmore. GOLDBAR: A Framework for Combinatorial Biological Design. ACS Synthetic Biology Article ASAP (2024). 

DOI: https://pubs.acs.org/doi/full/10.1021/acssynbio.4c00296

#
### Credits
The GOLDBAR-Generator was developed by [James Roberts]([https://github.com/Jamesr787]) at [CIDAR LAB](https://www.cidarlab.org) under [Douglas Densmore](https://www.cidarlab.org/doug-densmore).

#
