# ME700-Final-Project 

## Table of Contents

- [Skills Used in This Project](#skills)
- [Problem Overview](#Reproducibility)
- [Analysis](#Analysis)
- [Conda Environment, Installation, and Testing](#install)
- [How to Use The Codes](#htu)
- [References](#references)

## Skills Used in This Project <a name="skills"></a>
**Setting Up a Conda/Mamba Environment**  
At the start, I only had a basic understanding of Conda environments. However, throughout ME700, I learned how to create isolated environments, manage dependencies, resolve conflicts, and ensure reproducibility, making my workflow much more efficient and organized.

**GitHub and Modular Coding**  
I learned how to write more modular and maintainable code by organizing my project into logical functions and modules. This not only improved my code quality but also enhanced my ability to track changes.

**Automation with Bash Scripting:**
Managing multiple datasets and running numerous scripts manually is inefficient and error-prone. Through this project, I developed skills in automation, writing Bash scripts that allowed me to run multiple Python scripts and submit jobs to the computing cluster automatically.

**Matrix Structural Analysis (MSA):**  
This project marked my first time implementing Matrix Structural Analysis (MSA) from scratch. I learned the theoretical foundations, including assembling global stiffness matrices, applying boundary conditions, and solving for displacements and forces. This gave me a deep understanding of how MSA works and how it can be applied to various structural problems.

**Validation of Results:**
Given that my project focused on comparing two different analysis methods (MSA and Finite Element Analysis), I gained experience in systematically validating results.This experience reinforced the importance of verification and validation in any computational project.
  
### Technical Correctness:   


## Reproducibility

### Step 1: Copy files from GitHub to your desired destination (

### Step 2: How to setup a mamba/conda environment and install dependencies
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
### Step 3: Running codes   
Copy and paste 
```
```
Becarefull
### Step 4: Look at the results
   

## Code Structure

<p align="center">
<img src="https://github.com/sarajahedazad/ME700-Final-Project/blob/main/figures/final_project_flowchart.png" width="500">
</p>

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
## Analysis
**Matrix Structural Analysis**   

**Finite Element Analysis**

## Wrapping Up Results <a name="wrapup"></a>


## Things to be careful about when reproducing the results


## Results of the Project <a name="results"></a>


## References









