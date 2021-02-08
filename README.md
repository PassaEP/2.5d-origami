# 2.5d-origami

From https://github.com/PassaEP/CSV2SVG. 

Bootstrapping the DNA origami visualization process by generating simple block shapes with the Cadnano2.5 API and generating their forms as SVG files. 

## Dependencies 
Need cadnano2.5, pysvg-py3, numpy.. Follow installation instructions for cadnano2.5 at https://cadnano.readthedocs.io/en/master/installation.html#install-python-3-6. 
Otherwise: 

`pip install pysvg-py3`
`pip install numpy`

## Usage 

### Features (as of 2/6/21) 
- Annotated dimensions on 2D slices 
- Can generate XY, YZ, and XZ views
- Can generate as many pseudo 3D views as anyone would want, for any given angles (at your own peril). 
- Honeycomb lattices on 3D view as well

### In-Progress Features
- Need to render handles properly
- Switch to Click for command line usage 
