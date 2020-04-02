# Instructions to contribute data publicly available under DOI

This assumes that the data is in [Zenodo](www.zenodo.org).
It should work also for other sources, but has not been tested.

1. Collect the information in the first cell of the [Jyputer notebook](https://github.com/NMRLipids/NMRlipidsVIpolarizableFFs/blob/master/DataContribtionScripts/AddData.ipynb):

>DOI="**Give the DOI here**" \
>def_file  = "**Give the definition file used to calculate order parameters here**" \
> \
>user_information = """ \
>DOPE test \
> \
>#NMRLIPIDS BEGIN \
>@SIM \
>@SOFTWARE=**Give the name of used simulation program here, e.g., Gromacs, NAMD, etc.** \
>@FF=**Give the name of the force field, e.g., CHARMM36** \
>@FF_SOURCE=**Give the source of the force field, e.g., CHARMM-GUI** \
>@FF_DATE=**Give the date of the force field, e.g., when downloaded from CHARMM-GUI** \
>@TRJ=**Give the name of the trajectory in the Zenodo repository, e.g., trr or xtc file in Gromacs** \
>@TPR=**Give the name of the structure file in the Zenodo repository, e.g., tpr in Gromacs** \
> \
>#NMRLIPIDS END

Example can be found from [here]().

2. Deliver this information to us by commenting the blog or some other way,
or if you want to analyze yourself go to step 3.

## Running analysis for the data

3. Fork this GitHub repository to your own GitHub account and local computer.
4. Copy the information made in step 1 into the first shell of
the [Jyputer notebook](https://github.com/NMRLipids/NMRlipidsVIpolarizableFFs/blob/master/scripts/AddData.ipynb),
and define working directory at your local computer in the second cell.
5. Run the AddData.ipynb notebook to calculate the order parameters.
6. Push the directory (which should have name like ”10.5281:zenodo.3557459.0”) in the NMR-
lipidsVIpolarizableFFs/Data/Simulations directory to GitHub and make pull request to the master branch in the NMRlipids project.
