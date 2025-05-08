# ME700-Final-Project 

**Skills Used** 
Setting up a conda/mamba environment: Before I just had a shallow familarity with using conda env, I was not really that familiar wih how to set those envs and install dependencies, but throuout the course I learned how to use them.
Github: I learned how to write more modular codes and how to organize these codes better. 
Automation: when you are working with several dataset and pieces of code, you cannot really go through running them manually one by one. I learned how imoportant aoutomation
Matrix structural analysis: In this semester, it was my first time implementing matrix structural analysis. it was 


### Reproducibility


**A protocol that we can really really trust**

How to 
```
module load miniconda
mamba create --name me700-final
mamba activate me700-final
mamba install -c conda-forge matplotlib
mamba install -c conda-forge fenics=2019.1.0
mamba install -c conda-forge h5py
mamba install -c conda-forge scipy
pip install pandas
pip install pygmsh
pip install numpy
```
The following files should be in the same folder:
```
```
Then you shpuld open a terminal when you are in the folder as 
The following files should be in the same folder:
```
```
Then you should open a terminal 

**Code Structure**

<p align="center">
<img src="https://github.com/sarajahedazad/ME700-Final-Project/blob/main/readme%20figures/final_project_flowchart.png" width="500">
</p>

**Preprocess:**
### Before Analysis   
**Mesh Generation**   
Mesh generation relies on the routines defined in `archstruct_meshgen_functions.py`. To recreate the full dataset (of which only a handful of samples are used here), simply run the `archstruct_meshgen_main.py` script.

**Generating CSV Containing Parameters**  
You can specify the list of parameters to iterate over when running batch jobs in `archstruct_jobscsvgen.py`. For example, to perform analysis on three samples with cell keys 44444444, 55555555, and 66666666 across three mesh critical lengths (2, 1, and 0.1), you should define:
```
strkey_lst = ['44444444', '55555555', '66666666']
lcar_lst = [2, 1, 0.1]
```
You can see the rest of the parameters that should be defined in `archstruct_jobscsvgen.py`.   
Note: It should be noted that most of the parameters in this csv file won't be used in MSA, and parts of them won't be used in FEA either if there is no perturbation. Even if there is perturbation, not all of them will be used (it depends on the type of perturbation)
### Analysis
**Matrix Structural Analysis**   

**Finite Element Analysis**

### Wrapping Up Results


### Things to be careful about when reproducing the results


### Results of the Project


### References









