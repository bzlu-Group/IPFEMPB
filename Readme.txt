This program, IPFEM-PB, is developed based on the toolbox Parallel Hierarchical Grid (PHG), whch is currently under active development at State Key Laboratory of Scientific and Engineering Computing of Chinese Academy of Sciences. Please refer to http://lsec.cc.ac.cn/phg/ to get more details about the installation and usage of PHG. Afer PHG has been successfully installed, you need to change the content of 'PHG_MAKEFILE_INC' in Makefile to the installation path of PHG on your computer. Then you can compile and run our program on a parallel computer by running 'work.sh'.

You can modify 'work.sh' and 'example.options' to adjust some parameters of our program:

/*Basic Input*/
mpirun -np 32 ./IPFEMPB #denotes running IPFEMPB using 32 MPI processes,
-epsilon_solute 2.0 #denotes the relative permittivity of the solute molecule,
-epsilon_solvent 80.0 #denotes the relative permittivity of solvent,
-concentration 0.10 #denotes the salt concentration with unit mol/L,
-temp 298.15 #denotes the temperature with unit K,
-pqr XXX.pqr #denotes the PQR file,

/*Advanced Input*/
-box_length 10 #denotes the size of the calculation domain: [-Domain_lenth, Domain_lenth]^3 unit: 10^-10 m
-box_size 1.0 #denotes the size of the initial mesh with unit 10^-10 m
-PB_type LPBE #LPBE or NPBE denotes the type of PBE (LPBE or NPBE) to be solved
-refine0 12 #denotes the number of times the initial grid is uniformly refined, each time for dichotomous refinement, it is recommended to choose a multiple of 3
-dof_type P1 #P2, P3... denotes the type of finite elements space
-fn_mesh XXX.mesh #denotes using the given mesh as the initial grid
-xfem_ctol 1.0 #denotes the tolerance used to limit the relative curvature of the surface with respect to the size of neighboring elements
-ls_order 2 #denotes the order of polynomial interpolation for Gaussian molecular surface
-vtk #export vtk-file: XXX.vtk

We also provide a tool for visualization of the molecular surface potential:
Get XXX.vtk and run "python Get_Surface_Potential.py XXX.vtk", you can get the visualization of the molecular surface potential file XXX_Surface_Potential.vtk

