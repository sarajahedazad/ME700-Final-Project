# ME700-Final-Project 

### Table of Contents
- [Problem Definition]( #def )
- [ME700 Applied Skills](#skills)
- [Technical Correctness](#techcorrect)
- [Reproducibility](#Reproducibility)
- [Code Structure](#codestruct)
- [Results](#Results)
- [References](#References)
- [Other Things to Talk About!](#others)

## Problem Definition <a name="def"></a>
Architected materials are fascinating because their geometry—not just their composition—determines their mechanical behavior. In my research, I focus on strut-based architected materials and how they buckle when compressed.

But it's not just if they buckle—it's how they buckle: what shape do they take? How does that relate to their internal layout? That’s the deeper question I’m exploring.

I previously simulated these behaviors in FEniCS, but validation has been a challenge. So I decided to use the matrix structural analysis code I developed in ME700 to cross-validate my results and build a more reliable workflow.

The goal is to compare FEniCS simulations—using the nonlinear Simo-Reissner beam theory—with results from my matrix structural analysis code. I’ll compute the critical load, critical displacement, and buckling shape (mode), then compare them across both approaches.

Doing this manually for a few samples is easy. But as my dataset grows, it becomes tedious. That’s why I plan to build automation to run comparisons, extract results, and perform mesh refinement automatically.

This will serve as a modular, testable framework that not only helps with this class, but directly supports my PhD research moving forward.   

<p align="center">
<img src="figures/gif_C.gif" width="300">
</p>

## ME700 Applied Skills <a name="skills"></a>
**Setting Up a Conda/Mamba Environment:**   
At the start, I only had a basic understanding of Conda environments. However, throughout ME700, I learned how to create isolated environments, manage dependencies, resolve conflicts, and ensure reproducibility, making my workflow much more efficient and organized.

**GitHub and Modular Coding:**  
I learned how to write more modular and maintainable code by organizing my project into logical functions and modules. This not only improved my code quality but also enhanced my ability to track changes.

**Automation with Bash Scripting:**   
Managing multiple datasets and running numerous scripts manually is inefficient and error-prone. Through this project, I developed skills in automation, writing Bash scripts that allowed me to run multiple Python scripts and submit jobs to the computing cluster automatically.

**Matrix Structural Analysis (MSA):**  
This project marked my first time implementing Matrix Structural Analysis (MSA) from scratch. I learned the theoretical foundations, including assembling global stiffness matrices, applying boundary conditions, and solving for displacements and forces. This gave me a deep understanding of how MSA works and how it can be applied to various structural problems.

**Validation of Results:**   
Given that my project focused on comparing two different analysis methods (MSA and Finite Element Analysis), I gained experience in systematically validating results.This experience reinforced the importance of verification and validation in any computational project.
  
## Technical Correctness: <a name="techcorrect"></a>    
I believe my work is technically correct because the core objective of my project is to validate the results obtained from finite element analysis (FEA) using a second method, Matrix Structural Analysis (MSA). To ensure accuracy, I used MSA as an independent verification method to cross-check the results from FEA.

**Validation Approach:**   
*Selection of Structures:* I selected 10 symmetric or semi-symmetric structures for testing (5 from each). Semi-symmetric structures are those that maintain their shape when rotated by 180 degrees. These structures are specifically chosen because their eigenvalues can become negative, a behavior that is not guaranteed for asymmetric structures.

**Comparison Metrics:**   

*Critical Displacement:* I compared the critical displacement values between the two methods.  

*Critical Load:* I checked if the critical loads calculated by FEA and MSA matched.  

*Eigenmode Shapes:* I visually compared the eigenmode shapes between the two methods, ensuring they were consistent.  

**Further Improvements:**    
To enhance the accuracy and reliability of the validation process, I plan to:

- Compare Results with Theoritical Solutions: Implement a basic column analysis using both MSA and FEA, and compare the results with theoretical solutions for direct verification.

- Explore Mesh Sensitivity: Conduct a mesh refinement study to evaluate how the FEA results converge. (This was limited due to memory issues with MSA.)

- Include Asymmetric Structures: Extend the validation to asymmetric structures and develop a systematic method for comparing results between MSA and FEA for these cases.

- Automate Eigenmode Comparison: Develop an automated approach for comparing eigenmode shapes, such as using interpolation and normalization techniques, instead of relying on manual visual comparison.


## Reproducibility
### Step 0: Go to your desired destination on SCC and open a terminal

### Step 1: Clone files from GitHub
Copy and paste these in the terminal.
```
git clone https://github.com/sarajahedazad/ME700-Final-Project.git   
cd ME700-Final-Project
```
`Note:` Be careful that when running the codes, you are in the correct fodler where the codes are. `cd ME700-Final-Project` makes sure you are in the correct folder `ME700-Final-Project`.   
### Step 2: Setup a mamba/conda environment and install dependencies  
Cooy and paste this in the terminal.
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
**OR** You can alternatively use conda instead of mamba. (load conda instead of miniconda)
Once you have gone through this process, next time you need that environment you can simply run the following:
```
module load miniconda
mamba activate me700-final
```
### Step 3: Running the Codes   
Copy and paste the following commands in the terminal on SCC. You have to wait until the jobs are finished. You will recieve an email when the first job starts and when the last job finishes.
```
qsub -N myjob_meshgen job_meshgen.sh
qsub -N myjob_jobscsvgen -hold_jid myjob_meshgen job_jobscsvgen.sh
qsub -N myjobarr_msa -hold_jid myjob_jobscsvgen jobarr_msa.sh 
qsub -N myjobarr_fea -hold_jid myjob_jobscsvgen jobarr_fea.sh
qsub -N myjob_msavsfea -hold_jid "myjobarr*" job_msavsfea.sh

```
Note: to remove the email option look into files `job_meshgen.sh` and `job_msavsfea.sh` and remove the line of codes that start with 


When you are done with the environment, you can deactivate it by typing this in the terminal: `mamba deactivate`
### Step 4: Look at the results
You can look at the generated data for inputs, and the output results in the following directories:  
**Mesh**   
Mesh are generated using `gmsh` and stored in the designated directory. 
*Directory:* `inputs/archstruct/mesh/mesh2D/meshlcar1`    (`1` is `lcar` or carachteristic length of the mesh here. We define it when we generate a mesh.)
*File Naming Example:*  `meshny4nx2celltypes456cellsize20lcar1ID1.xdmf` (Each `.xdmf` file comes with a `.h5` file. )

**Keys**
*Directory:* 
*Naming Example:* `keys_ny4nx2celltypes456.csv` 

**MSA Results: Critical Displacement, Critical Force and Eigenvector:**   
*Directory:* `outputs/archstruct/results2D/FEA/resultslcar1`     
*Naming Example:* ``, ``, ``    
Note: The eigenvectors that are saved for MSA results, include the terms that are related rotation, unlike the eigenvector array that is saved for FEA that only include terms related to x and y displacements. The eigenvector saved for MSA is 1D, while it is 2D for FEA.    

**FEA Results: Critical Displacement, Critical Force and Eigenvector:**
*Directory:*    
*Naming Example:*    
Note: The eigenvectors that are saved for MSA results, include the terms that are related rotation, unlike the eigenvector array that is saved for FEA that only include terms related to x and y displacements. The eigenvector saved for MSA is 1D, while it is 2D for FEA.      
**Wrapped Up Results:**


## Code Structure  <a name="codestruct"></a>

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
### Analysis
**Matrix Structural Analysis**   

**Finite Element Analysis**

### Wrapping Up Results <a name="wrapup"></a>


### Things to be careful about when reproducing the results


## Results  
**Comparison of Critical Displacemnt Values in MSA and FEA Across Ten Samples**   
<p align="center">
<img src="https://github.com/sarajahedazad/ME700-Final-Project/blob/main/figures/jobssamplecheck_dispcritical.png" width="700">
</p>

**Comparison of Critical Force Values in MSA and FEA Across Ten Samples**   
<p align="center">
<img src="https://github.com/sarajahedazad/ME700-Final-Project/blob/main/figures/jobssamplecheck_forcecritical.png" width="700">
</p>

**Percentage Errors**   
<p align="center">
<img src="https://github.com/sarajahedazad/ME700-Final-Project/blob/main/figures/jobssamplecheck_errors.png" width="700">
</p>


**First Eigenmode Comparison**  

It should be noted that the configuration for each analyasis was saved individually and then put together in one figure using illustrator.  
<p align="center">
<img src="https://github.com/sarajahedazad/ME700-Final-Project/blob/main/figures/msavsfea_eigenmodes.png" >
</p>

**Conclusion**  
The percentage error for critical loads and displacements across 10 samples ranges from approximately 0.5% to 2.3%. The eigenmode visuals are also nearly identical. These results could potentially improve with mesh refinement.

## References  
* ChatGPT   
* ME700 Course Material   
* [Arclength Displacement Control Tutorial](https://github.com/pprachas/fenics_arclength/blob/master/examples/displacement_control/fiber_network.ipynb)   


## Other Things to Talk About!  <a name="others"></a>
* If you really hate having the job logs in your folder (like me!) just type `rm -f myjob*` in your terminal after you are done with getting the results. 









