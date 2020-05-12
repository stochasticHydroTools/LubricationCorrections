# LubricationCorrections
Fast algorithms for applying lubrication corrections to fluctuating suspensions of spherical particles above a bottom wall

In order to run the code, the SuiteSparse library must be installed on the system
`http://faculty.cse.tamu.edu/davis/suitesparse.html`

In Ubuntu 18, this can be accomplished by running
`sudo apt install libsuitesparse-dev`

Additionally, the C++ libraries `Boost` and `Eigen` must be installed
Ensure that the verion of boost is 


In addition, the following python packages must be installed via pip

```
pip install --user scipy
pip install --user pyamg
pip install --user scikit-sparse
```
To run the code, one needs to `Make` the C++ helper code in the directory `/Lubrication`


```
https://github.com/stochasticHydroTools/LubricationCorrections/tree/master/Lubrication
```
