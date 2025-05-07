# ME700-Final-Project

I am work shppoing my ideas here!  
I can use this project as an opportunity to really organize my own reseach.

**What do I wanna add to my project?** 
- simple column vs theory --> different meshes
- 2D matrix structural analysis vs fenics -> different meshes
  or now do it for symmetric structures. Do everything for symmetric ones.
if I have the time, I have to add the asymmetric ones

***Documentation***
- Theory  
Maybe I can take this time to understand the concepts better? (this might take like 6 hours)
(2 hours for theory, 2 hour sfor reviewing fem and all, 2 hours for the code)
- Actually learn the kind of chnages to MSA (2-3)
- How to run codes documentation (2 hrs)   

- Mesh generation clenaning (2 hr)
- Functions cleaning (2 hr) 
- Automation tools in one place (2 hr)

6-3-8  

**Skills Used**   
ADD THE SKILLS YOU USED IN THIS PROJECT HERE

**Reproducibility**   
ADD REPRODUCIBILITY   

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









