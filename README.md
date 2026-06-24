# MESA-DP
MESA code to simulate the impact of dark photons on stars

This code is used in the publication arXiv:2606:xxxxx. The physics and technical details are explained in this paper. Please cite it if you use the code in your research.

# MESA Installation
To use the code, you first need to install MESA and MESA SDK following the instructions here:

[1] https://docs.mesastar.org/en/r15140/installation.html

Some useful tutorials: 

[2] https://mesa-leuven.4d-star.org/tutorials/tuesday/morning-session/

[3] https://jschwab.github.io/mesa-2021/solutions.html#org417d2fb

[4] https://jschwab.github.io/mesa-2017/solutions.html#orgf3d2cdc

## A simple guide to the installation
 - Download `mesa-25.12.1.zip` and `mesasdk-x86_64-linux-25.12.1.tar.gz` (version numbers may change) from [1]; unpack the two packages anywhere. For example:
  ```bash
  cd .../path_to_MESA
  unzip mesa-25.12.1.zip
  tar -xzfv mesasdk-x86_64-linux-25.12.1.tar.gz
  ```

 - Modify ~/.bashrc or ~/.bash_profile to add the following lines:
   ```bash
    # set MESA_DIR to be the directory to which you downloaded MESA
    # The directory shown is only an example and must be modified for your particular system.
    export MESA_DIR=.../path_to_MESA/mesa-25.12.1
    
    # set OMP_NUM_THREADS to be the number of cores on your machine
    export OMP_NUM_THREADS=2
  
    # you should have done this when you set up the MESA SDK
    # The directory shown is only an example and must be modified for your particular system.
    export MESASDK_ROOT=.../path_to_MESA/mesasdk
    source $MESASDK_ROOT/bin/mesasdk_init.sh
    
    # add shmesa (the MESA command line tool) to your PATH
    export PATH=$PATH:$MESA_DIR/scripts/shmesa
   ```
 - Restart the terminal (or open a new terminal), run the following commands (no sudo is required):
   ```bash
      cd $MESA_DIR
     ./install
    ```
 - If  mesasdk_init.sh complains "missing prerequisites:  C shell", then you may need (on Debian/Ubuntu)
    ```bash
    sudo apt install tcsh
    ```
# How to use this code
 - Create a MESA project using from the template:
  ```bash
  cp -r $MESA_DIR/star/work mesa-DP
  cd mesa-DP
  ```
 - Copy the files downloaded from this git repository into mesa-DP

 - Compile the new files
  ```bash
  cd mesa-DP
  ./mk
  ```
 - If you just need to simulate one star with a specific profile (including new pyhsics settings), simply run 
  ```bash
  cd mesa-DP
  ./mk
  ```
 - If you need to run the simulation for many stars, please have a look at `run_models_v2.1.py`.
 - If you need to run the simulation for the Sun with many sets of new physics parameters, please have a look at `run_models_Sun.py`.



