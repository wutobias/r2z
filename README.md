## Z-matrix conversion with RDkit molecules

With this little python implementation, RDkit molecules can be converted to a Z-Matrix topology representation. If the RDkit molecules holds conformations, one can also generate the Z-Matrix coordinates or convert the Z-Matrix coordinates back to Cartesian coordinates.
During all conversions one easily can make a mistake by using wrong units, there the `openmm` is used for dealing with units and enforcing the corrects during every step.

## Coordinate transformations
The transformation from Cartesian coordinates is straightforward. However, the back conversion is more tricky and somewhat error prone. One difficulty in this context is to avoid numerical instabilities, for instance through round-offs, since these might propagate through the Z-Matrix during the conversion. In order to circumvent this, we implemented the Natural Extension Reference Frame algorithm, which minimizes these errors. The algorithm is described in detail in these two articles https://doi.org/10.1002/jcc.20237 and https://doi.org/10.1002/jcc.25772

A general remark on transformations is that the conversion from ZMatrix coordinates to Cartesian coordinates requires 3 reference Cartesian coordinates (e.g. coordinates in a host or protein molecule) and 3 torsion angles, 2 angles and 1 bond length with respect to these coordinates. These coordinates are also called 'virtual coordinates'. For instance, let's say we have the Zmatrix coordinates for a molecule (e.g. a ligand) with 5 atoms: `A-B-C-D-E` and we want to convert them back to Cartesian space. First, we need 3 Cartesian coordinates `X-Y-Z` in the reference frame, which can be basically *anything* in the lab coordinate system, as long as it is well defined. Second, we need the torsion angles `X-Y-Z-A`, `Y-Z-A-B`, `Z-A-B-C`, the angles `Y-Z-A`, `Z-A-B` and the bond length `Z-A`. By convention, the Cartesian coordinates of our molecule `A-B-C-D-E` will be shifted to the position of `Z`. The calculation of all these virtual coordinates can be carried out with this python class.

## Requirements
* openmm
* rdkit
* numpy

## Examples
Different examples for how to use this python class, can be found in examples.

## Remark
The algorithm is quite inefficient for large molecules. Peptides are still ok, but computing proteins might take a while.
If the RDkit molecule object contains more than one molecule, you won't be able to build the ZMatrix. In case you want to process more than one molecule (e.g. a host guest complex), either add a bond between them or build the Z matrices for each molecule separately. The second option is the cleanest procedure and would require to build a reference frame from the first molecule that is used for transformations of the second one.
