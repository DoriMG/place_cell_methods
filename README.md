# Place cell methods
Code accompanying manuscript on performance of place cell methods. This code allows the user to create datasets of model place cells. Several features of the place cell can be altered:
* the place field peak
* the place field width
* variability across traversals
* reliability across traversals
* noise in the signal

The default version imports a cell array (`allTras.mat`), with each cell containing the location over time in a single traversal, obtained from real mice running through a virtual maze. Users can opt to import their own location data in this same format to be used to create the model dataset.

## Usage
Create a model dataset with the default settings using:

`data = model_place_cells()`
![Rasterplot of example data](/images/example_dataset.png)

Create a model dataset with variable place cells with a place field width of 100 cm:

`data = model_place_cells('bin_size', 20, 'var_sigma', 0.5)`
