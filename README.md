 
# Mani Research Group

Mechanical Engineering Dept., Stanford University


## Channel Flow Code

**Date created : 05/07/2008**
by Jongmin Seo & Sanjeeb T. Bose 

**Date modified : 04/23/2018**
by Danah Park


---

### Instruction to the code
**Date writen : 01/24/2019**
by Danah Park


0. Intro-introduction
   * This script is generated to help people who use this code. Understanding this code can be very time-consuming and stressful (at least it was for me). Hope this script is helpful! :)


1. Introduction
   * This is a computational fluid dynamics solver for a 3D turbulence channel.
   * The Inverse Macroscopic Forcing Method is implemented to the code.


2. Code compilation and monitoring
   * The code uses MPI Fortran for the compilation.
   * Libraries it uses are MPI, MKL, LAPACK, FFT by Sandia etc, (see src/makefile)
   * Step by step compilation
       * Compile the code at `src/` using `make`. This will make a binary file `bin/channnel`.
       * Set input parameter in `input/input_io.dat`.
       * Change slurm settings in `tools/DNS_slurm.sh`. Remeber to set `WORKDIR` to your working directory.
       * At `tools` run `bash run_dns`. This will submit your job.
   * Code output and monitoring
       * When you run the code there will be output files, job monitoring log file, and `diagnostics.dat` file.
       * Output files are stored in `output/flow_field_*/` in binary data. See 4 for the data structure of each data file.
       * Log files are stored in your working directory. `log*.out` file will tell you how code is running and `log*.err` will tell you if there is any error in the code.
       * `diagnostics.dat` is the most important file to detect running errors. In is in the order of time step, eddy-turnover time, CFL number, u_\tau for the bottom wall, u_\tau for the top wall, v_\tau for the bottom wall, v_\tau for the top wall, transporter forcing term at a certain location near wall.
           * CFL number should be less than 1
           * u_\tau for the bottom wall, u_\tau for the top wall is around 1 when statistically stationary
   * Code and job cleaning
       * You may monitor your job with `squeue` and kill your job with `scancel 'jobid'`.
       * Use `make clean` to delete compile objects and execution file
  
  
3. Program directory and modules
    * Program directories
        * `input/`: stores input data
        * `output/`: stores output data
            * `output/flow_field`: stores transporter velocity (u) and pressure
                * Each snapshots are stored at eddy-turnover time of the file name. For example, `field_101_700.dat` file is a snapshot of the flow field at eddy-turnover time of 101.700. Same goes for below two.
            * `output/flow_field_transportee`: stores transportee velocity (v) and pressure
            * `output/flow_field_IMFM`: stores transporter forcing term
            * `output/mean_profile`: stores transporter forcing term
            * `output/restart`: stores restart file
                * You can set restart file using flow_field snapshots. Just change names of field data in to `.in` file and store it in this location. For example, `output/flow_field/field_101_700.dat` to `output/restart/restart.in`, `output/flow_field_transportee/field_101_700.dat` to `output/restart/restart_transportee.in` and `output/flow_field_IMFM/field_101_700.dat` to `output/restart/restart.in`. If you are not using IMFM, `restart_IMFM.in` is enought to have to restart your computation.
        * `src/`: stores code
        * `stats/`: stores postprocessing modules
        * `tools/`: stores Slurm scripts and bash scripts
        * `bin/`: stores the binary execution file
    * Program modules in `src/`
        * `comm_routines.f90`: generates communicator of MPI
        * `fftw_f77.i`: stores FFTW parameters
        * `global.f90`: controls global variables, including allocation and deallocation
        * `grid.f90`: generates grid
        * `main.f90`: main function that call other modules
        * `makefile`: makefile
        * `numerics.f90`: numerics for u are all done here including discretization and time advancing
            * I only used the code with AB/CN as time discretization and CD2. The code and MFM implementation may or may not run on other schemes.
        * `numerics_tee.f90`: numerics for v are all done here including discretization and time advancing
        * `output.f90`: read input files, generates output files, dianostics file and etc.


4. Data structure
    * Output data structure for `output/flow_field` and `output/flow_field_transportee`
        * Output data is stored in the order of u1, u2, u3, and p in Fortran order (this is VERY important! Fortran uses column major order unlike C, python and all other languages except fortran). I attached my snapshot reading python code here. This code will read a snapshot and change it into u1, u2, u3 in a row major order.
        
        ```python
        # Data structure of the snapshot file
        dt = np.dtype([('id', np.int32), ('time', np.float64),
                       ('u1', np.float64, (Nx*Ny*Nz)), # note that this is a double type
                       ('u2', np.float64, (Nx*Ny*Nz)),
                       ('u3', np.float64, (Nx*Ny*Nz)),
                       ('p', np.float64, (Nx*Ny*Nz)),
                      ]) 

        # Read the snapshot in binary file and change it something that we can read
        a_f = np.fromfile('field_101_700.dat', dtype=dt)

        # Save the data read
        t=a_f['time'][0]
        u1=np.array(a_f['u1'][0])
        u2=np.array(a_f['u2'][0])
        u3=np.array(a_f['u3'][0])

        # Reshape into 3-dimensional array
        u1 = np.reshape(u1, (Nx, Ny, Nz), order='F')
        u2 = np.reshape(u2, (Nx, Ny, Nz), order='F')
        u3 = np.reshape(u3, (Nx, Ny, Nz), order='F')
        ```


5. Numerics
    * See Jongmin's thesis... I should go get lunch..


6. Postproccessing
    * The postprocessing codes are developed using python and stored in `stats/`. How each modules work are explained at the beginning of each .py file. You can develop MPI Fortran or C++ code to do the post-processing or use higher level code like Matlab or python to postprocess. Since the output files are all binary, you need to have some sort of post-processing code just to read one snapshot data. If you are using what is already developed please check when using stat modules other than `uv_calculation/`, `u_calculation/`, and `v_calculation/`. These three modules are well-verified and others are not as so.


7. Caveats
  * This code was first built as DNS, LES solver for channel and later was implemented to superhydrophobic application. But the solver is no longer able to do either, only has some remnants of the previous implementations. Later, I have implemented the Macroscopic Forcing Method with transporter variables. Hence, the code is hard to read.
  * When you are using Certainty cluster partition `all` or `24cores` to run your job it needs exactly '24 mutiples' of threads to run. Lower level issue...
  * `output/mean_profile` just ignore the files in here. It gives wrong data.
