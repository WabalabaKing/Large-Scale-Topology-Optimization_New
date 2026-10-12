# Short Cantilever Beam (SCB) Beam Problem
## Problem Description 
This is a demonstration of Compliance minimization with volume fraction constraint and CG  constraint.
## File Structure
The files are structured as follows:
```
SCB
|-Mesh
|  |-SCB.nam
|  |-RunFilt.inp
|-SCB_CG_5
|  |-RunOptim.py
|  |-MBB.inp
|-SCB_CG_9
|  |-RunOptim.py
|  |-MBB.inp
```
Each `SCB_*` is the root directory of the respective optimization problem.
`Mesh` is the location where the mesh and filtering matrix will be stored.

## Running The Problem
### Filtering
The first step is to construct the filtering matrix for the topology optimization problems. 
1. Navigate to  `Mesh`
2. Open `RunFilt.inp` and make sure that the second line points to the location of the mesh file.
3. Run `calFilt.exe <INPNAME> -r <RMIN> -f <NNZ> -Passive <0 or 1>`
   
   3.1 For this case, run `calFilt.exe RunFilt -r 0.04 -f 1000 -Passive 1`
   
   3.2 For parallel filtering, run `export OMP_NUM_THREADS=<NUMofTHREADS>` prior to running `calFilt`.
   
   3.3 Once done, `dnnz.bin, drow.bin, dval.bin, dcol.bin, dsun.bin` will appear in the mesh directory.
   
   3.4 The argument `-Passive` argument will govern if the "passive element will be excluded during filtering". If `-Passive 1` ONLY active elements are considered during filtering, i.e. there will be a clean cut between the skin and internal structure. If `-Passive 0` then passive elements are considered during filtering, there will be a layer of partial density elements at the passive-active element interface due to filtering.
   

### RunOptim.py
1. Before running the problem configure `RunOptim.py`. The input section is between the `Line22-Line41`. A few things to note is MeshDir should be the ABSOLUTE directory pointing to the Mesh directory that we are working with in the previous section.
The `objectives.csv` is formatted as such:
  `COMPLIANCE	 FULL SOLID VOLUME	 DESIGN VOLUME	 VOLUME_FRACTION	 DISCRETENESS	 MASS	 CGx	 CGy	 CGz	 PNORM` 
  And therefore `TableReader(0,3...)` points to `VOLUME_FRACTION` and `TableReader(0,0...)` points to `COMPLIANCE`
2. The optimizer setup is toward the later part of the code using `nlp.set()`. For more references, checkout https://coin-or.github.io/Ipopt/OPTIONS.html

### SCB.inp 
  The `*.inp` file are written similar in style of NASTRAN input files. Before running, make sure the 2nd line `*INCLUDE, INPUT=` points to the absolute location of the mesh file in `\Mesh`

### Running:
  After modifications and configurations outlined by the previous steps, run:
  `python RunOptim.py >log.opt &` to run the optimization problem. The `>log.opt` is optional, but very helpful in observing optimizer behaviour during the run.
  Similarly, you may run `export OMP_NUM_THREADS=<NUMofTHREADS>` before to use multiple threads.

## Results:
After running both problem, another folder `FinalCounter1` should appear. (If you are doing continuation, a `FinalCounter*` should appear for each step taken during continuation. This directory stores the last iteration's output from `calTop.exe`. Inside this director `FinalCounter1\Direct` There are the following that is worth investigating:
```
elastic_Field_active.vtu is the post-processing of the design domain
objectives.csv is the summary of objective value and constraint values
density.dat is the material distribution. It might be useful when doing other forms of post-processing or restarting the problem. It has to be used together with the meshes and skinelements as it constains no geometry definition and element index information
```

The two cases with different CG constraints solved using LOQO are shown below. It was observed that for with a xCG constraint of 9m, the material distributrion is further away from the fixed edge.

<img width="1433" height="508" alt="image" src="https://github.com/user-attachments/assets/0f1dc499-5595-4d64-b4e0-a5dc702eddfc" />
