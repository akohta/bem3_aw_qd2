# bem3_aw_qd2

This is the three-dimensional acoustic wave scattering analysis program for two-dimensional periodic arrangement objects irradiated by a plane wave. 
The radiation force acting on the object can be analyzed. 
This is based on boundary element method, the own developed numerical solution is used. The corner problem free.
Intel Math Kernel Library (MKL), Gnu scientific library (GSL) and libpng are required. 
Gmsh is used to create mesh data for objects. 
The calculation program of quasi-periodic Green's function "d3_qpgf_d2" is used. 

![analysis model](model_qpbc2.png "analysis model (model_qpbc2.png)")  

## Usage of example code  

1. type 'make' command to compile.  
   The executable aw_d3qd2_bv_solver, example1.out, example2.out, example3.out are created. 
   The aw_d3qd2_bv_solver is the main solver of boundary integral equations. 
   The example1.out is the executable of source code example1.c, it shows a simplest example using "bem3_aw_qd2". 
   The example2.out is the executable of source code example2.c, it shows a example of sound pressure intensity analysis. 
   The example3.out is the executable of source code example3.c, it shows a example of outputting the instantaneous value of sound pressure as an image.  
   
2. type './aw_d3qd2_bv_solver' with arguments of plane wave datafile name, periodicity datafile name, medium datafile name, mesh datafile name, output datafile name, rotation and translation settings (optional).   
   For example, './aw_d3qd2_bv_solver pw.txt periodicity.txt medium.txt sphere_m1.msh ex.dat'. 
   The pw.txt is the sample of incident field datafile, a plane wave is defined in it. 
   The periodicity.txt is the sample of periodicity datafile, periodic boundary condition and lattice vectors are defined in it. 
   The medium.txt is the sample of medium datafile, one medium is defined in it. The domain number is assinged to the medium from 1 in order. 
   The sphere_m1.msh is the sample of mesh datafile, it is a sphere object. 
   It was created by using Gmsh geometry file sphere_m1.geo in the mesh_sample folder.
   The sphere_m1_image.png is the visualization result of the sphere_m1.msh. 

3. type './example1.out' with an argument of datafile name output by aw_d3qd2_bv_solver.  
   For example, './example1.out ex.dat'. 
   This executable calculates sound pressure, particle velocity, radiaton force and torque.  

4. type './example2.out' with an argument of datafile name output by aw_d3qd2_bv_solver.  
   For example, './example2.out ex.dat'. 
   This executable calculates sound pressure intensity distributions, outputs them to text files. 
   The I_example2.png is the visualization result of the intensity distributions, created by using Gnuplot script gscritp_example2.plt
   (using ImageMagick to convert eps to png).  
   
5. type './example3.out' with an argument of datafile name output by aw_d3b2_bv_solver.  
   For example, './example3.out ex.dat'. 
   This executable calculates instantaneous value of the sound pressure, outputs them to png image files. 
   The image files are output to the folder which has a name adding "images" to the datafile name specified in the argument (file-extension is excluded). 
   Each image file has a name that indicates the cross section, field component and number of time steps (ex. xz_p_014.png). 
   The color bar is output as color_bar.png in the same folder. 
   The range of color bar in each cross section is output to the info.txt file (ex. xy_info.txt for z=0 plane). 
   The xz_p.gif, yz_p.gif and xy_p.gif are animated gifs that concatenate the png files created by using the shell script gif_animation.sh.  
   
Please see d3qd2_src/bem3_aw_qd2.h for detail of functions. 
The main parts of the code are parallelized by using OpenMP. 
The number of threads is controlled by the environment variable OMP_NUM_THREADS.  

![mash image 0](sphere_m1_image.png "mesh image of the object (sphere_m1_image.png)") 
![intensity distributions 0](I_example2.png "intensity distributions (I_example2.png)")  
![xz_p.gif](xz_p.gif "instantaneous value of the p on y=0 plane (xz_p.gif)")![yz_p.gif](yz_p.gif "instantaneous value of the p on x=0 plane (yz_p.gif)")  
![xy_p.gif](xy_p.gif "instantaneous value of the p on z=0 plane (xy_p.gif)")  


## Analysis sample 1 (in the analysis_sample1 folder)  

This is the analysis result when arranged in a triangular grid pattern. 
The others are the same as the above example.

![intensity distributions 1](analysis_sample1/I_example2.png "intensity distributions (analysis_sample1/I_example2.png)")  
![xz_p.gif 1](analysis_sample1/xz_p.gif "instantaneous value of the p on y=0 plane (analysis_sample1/xz_p.gif)")![yz_p.gif 1](analysis_sample1/yz_p.gif "instantaneous value of the p on x=0 plane (analysis_sample1/yz_p.gif)")  
![xy_p.gif 1](analysis_sample1/xy_p.gif "instantaneous value of the p on z=0 plane (analysis_sample1/xy_p.gif)")  


## Verification  

The verification results are in the folder verification. 
The analysis result of symmetrically arranged 25 spheres using "aw_msp_ivf" (in the folder aw_msp_ivf_result) 
and the analysis result of periodic arrangement of spheres are shown. 


## About mesh file

This code can use quadrangular (bi-linear) and triangular (linear triangular) elements. 
I recommend using quadrangular element for reduce required memory. 
The samples of mesh data are in the folder mesh_sample. 
The file with extension .geo is the Gmsh geometry file. 
The file with extension .msh is the mesh file created by using Gmsh geometry file. 
These mesh files are created by the command 'gmsh -2 -tol 1.0e-15 xxxx.geo' in command line (xxxx.geo is a geometry file). 
The domain number (Physical Surface) 99 is assigned to the open region in Gmsh geometry file, because Gmsh can't use the number 0 (assigned to open region in the code). 
Please refer to the manual of Gmsh for detail of geometry file.


## Reference

1. Intel Math Kernel Library [MKL](https://software.intel.com/mkl)  
2. GNU Scientific Library [GSL](https://www.gnu.org/software/gsl/)  
3. The official PNG reference library [libpng](http://www.libpng.org/pub/png/libpng.html)  
4. Three-dimensional mesh generator [Gmsh](https://gmsh.info/)
5. The command-line driven graphing utility [gnuplot](http://www.gnuplot.info/)  
6. The utilities for manipulating images [ImageMagick](https://imagemagick.org/)  
7. The calculation program of quasi-periodic Green's function [d3_qpgf_d2](https://github.com/akohta/d3_qpgf_d2/)
8. The acoustic wave scattering analysis program [aw_msp_ivf](https://github.com/akohta/aw_msp_ivf)  
