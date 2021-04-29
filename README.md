# Registration
Spatial and temporal 3D image registration for mesoscale studies of brain development

USAGE:

1. Add the folder with the scripts to the Matlab path.

To do so, open Matlab and navigate in Matlab to where the folder `Registration` is stored. Right-click on the folder `Registration` and, from the dropdown menu, select 
`Add to Path`, then `Selected Folders and Subfolders`.

2. Import the brain images

In Matlab, navigate to a folder containing the z-stack images of a single 3D brain sample. Once in the folder, execute the following command in the Matlab Command Window:

```matlab
IM = mload('*.tif', 3, 1/4);
```
The first parameter here is a regular expression, which allows to select the files from the current folder. Here `'*.tif'` means that we select files with any name (`*`) that are in the TIFF format (`.tif`). We select all files because the folder we are currently in contains only the files of a single brain sample. More information on regular expressions can be found at https://en.wikipedia.org/wiki/Regular_expression

The second parameter here (`3`) indicates how much we coarsen the image resolution, compared to the original resolution of the TIFF images. For example, if the original XY-resolution in TIFF images was equal to 4 um/voxel, the final resolution in the images loaded into Matlab would be 4 * 3 = 12 um/voxel, which is an optimal resolution for our algorithm. We coarsen the resolution to average out individual cells in the images and to deal with brain areas instead. Once registration is finished, its results may be applied to the full-resolution images, if necessary.

The last parameter, `1/4`, shows the ratio between the XY-resolution within the original TIFF images and their Z-resolution, i.e. the distance between the neighboring TIFF images. For example, if the XY-resolution is equal to 4 um/voxel and the Z-resolution is 12 um/voxel, then the ratio between them would be equal to 4/12 = 1/4, as shown in the example. This parameter is needed to load the 3D images of the brain with the same XY- and Z-resolutions, defined by the previous parameter (e.g. 12 um/voxel in each direction).

The output parameter `IM` defines the name, under which you would like to load the brain sample. Choose unique names for different brain samples to avoid overwriting. It is recommended to save the brain images after each operation. To do so, once all brain samples are loaded, select them in the Workspace panel in Matlab, right-click on them and select the `Save` option from the dropdown menu.

3. Preprocess the loaded brain images

In Matlab, navigate to the folder `Registration` and run the following command:

```matlab
IMs_preproc = mprepare({IM1, IM2, IM3});
```

Here the input parameters in curly brackets (`{IM1, IM2, IM2}`) indicate the names of the brain samples to be preprocessed. It is important to preprocess all brain samples at once so that the resulting preprocessed samples (written to the new variable `IMs_preproc`) share the same size, orientation, and brightness settings.

The script `mprepare` is interactive. For each brain sample, you would be required to go through two steps. First, a window will pop up, containing four XY-slices of the 3D brain at different depths, and a color bar on top. By clicking on different areas of the color bar, choose the image intensity threshold such that most of the image background has zero intensity (color-coded blue) but the entire brain sample retains non-zero intensity (color-coded red/yellow). Having a zero-intensity background is important for the next steps. Once the optimal threshold is selected, click anywhere outside the color bar to save the choice.

Shortly after, another window will pop up, asking to choose the orientation for the brain. For the first brain sample, you may choose any orientation you prefer; for the following brain samples, the orientation should match that of the previous sample. The current brain is shown in green, whereas the previous brain is shown in red. To choose the orientation, click on any of the three projections of the brain sample in the image. The projection you selected will be flipped in both directions and the updated image will be displayed. Once the orientation is selected, click anywhere outside the image to save the choice. Save the preprocessed images `IMs_preproc`, similarly to how you saved the imported brain images on the previous step.

4. Make the brain samples symmetrical (optional)

To make the images of the brain samples symmetrical with respect to the center plane of the 3D image, run the following command:

```matlab
IMs_symm = mselfalign(IMs_preproc, 'young');
```
Here `IMs_preproc` is the variable containing all images of the brain samples preprocessed on the previous step; `'young'` is an age flag used to apply the optimal settings for the `'young'` or the `'adult'` brains; `IMs_symm` is the output variable containing the symmetrized brain images. The expected algorithm's runtime is around 5 minutes per brain sample, during which the progress will be displayed.

5. Register the brains to common coordinates

To register the brains, run the following command:

```matlab
[IMs_registered, AGES_ADJ] = mwarp(IMs_symm, [1 1 2], 'young', 'full');
```
Here `IMs_symm` is the variable containing all images of the brain samples symmetrized on the previous step (alternatively can be `IM_preproc` if the previous step was omitted). `[1 1 2]` is an array of the ages of the brain samples (i.e. the first brain loaded on `Step 3`, namely `IM1`, is a P1 brain; `IM2` is also a P1 brain; `IM3` is a P2 brain). `'young'` is an age flag used to apply the optimal settings for the `'young'` or the `'adult'` brains. `'full'` is a size flag indicating whether to align the `'full'` brains or their `'halves'` (hemispheres).

Aligning the `'halves'` requires symmetrizing the brains on the previous `Step 4` and allows to double the size of your dataset by superimposing left and right hemispheres. `IMs_registered` is the output variable containing the brain images registered to common coordinates. `AGES_ADJ` are the adjusted (developmental) ages of the brain samples in `IMs_registered` computed using the similarities of the brain samples (e.g. `[0.8, 1.1, 1.9]`). The expected algorithm's runtime is around 5 minutes per brain sample in the case of the `'full'` brains and around 3 minutes per brain for the `'halves'`.

6. Display the dynamics of brain development.

To evaluate the dynamics of brain development, run the following command:

```matlab
mmorph(IMs_registered, AGES_ADJ, 'full');
```
The command will reconstruct the dynamics of brain development, showing how the image intensity changes over time. In a separate window, it will also highlight the regions where the image intensity has increased (color-coded red) or decreased (blue) over the course of the last day. To save these results, the command will produce two videos, `growth_mean.avi` and `growth_diff.avi` respectively. Besides, it will produce a set of images with similar names, showing the average image intensity and its changes for each day within the range of the adjusted (developmental) ages `AGES_ADJ`.

7. Export the registered brain images (optional)

To export a registered brain image as a series of TIFF images, create a new folder, enter it in Matlab and run the following command:

```matlab
msave(IMs_registered{1});
```
Here `IMs_registered{1}` is the first of the registered images `IMs_registered`, i.e. `IM1` loaded on `Step 3`.
