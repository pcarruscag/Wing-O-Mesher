# Wing-O-Mesher
A half-plane mesh generator / morpher for wings (O-Grid closed into hemisphere).

## Description
The wing geometry is obtained by extruding a 4 digit profile and then applying linear distributions of taper, sweep, and twist.
Mesh morphing (using RBF's) can then be applied to obtain more complex geometries.
The mesh is a 2D O-Grid extruded along the span and then closed such that the farfield forms an hemisphere.
The wing itself can also be meshed to create meshes for FSI problems.
Meshes (2-D or 3-D) are exported in SU2 native format (ASCII).

This code is provided for transparency to simplify the eventual reproduction of some published work.

## Compilation
RBF features are implemented (efficiently!) in C++, they are setup to compile with [CodeBlocks](https://www.codeblocks.org/).
The rest is implemented in Octave/Matlab.

## Usage
`cp example_settings.txt mesh_settings.m` modify as needed and `octave RUNME.m`.

## License
[GPL-3.0](https://www.gnu.org/licenses/gpl-3.0.html)

