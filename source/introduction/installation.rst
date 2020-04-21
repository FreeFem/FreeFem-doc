.. role:: freefem(code)
   :language: freefem

.. role:: bash(code)
   :language: bash


Installation guide
==================

To use FreeFEM, two installation methods are available: user access (binary package) and access developers (from the source code).
Follow the section corresponding to your type of installation.

.. note:: Since the version 4.5, FreeFEM relese provides with the last version of PETSc.

|
|

Using binary package
--------------------

First, open the following web page :ref:`download page <download>` and choose your platform: Linux, MacOS or Windows.

.. note:: Binary packages are available for Microsoft Windows, MacOS and some Linux distributions. Since the release 4.5, FreeFEM binaries provide with the current version of PETSc.

Install **FreeFEM** by double-clicking on the appropriate file. Under Linux and MacOS the install directory is one of the following ``/usr/local/bin``, ``/usr/local/share/freefem++``, ``/usr/local/lib/ff++``

Windows installation
^^^^^^^^^^^^^^^^^^^^

.. note:: The windows package is build for Window 7 64bits. The support ended for all releases under Windows 32 bits since the V4. 

First download the windows installation executable, then double click to install **FreeFEM**.  
Install MSMPI for parallel version under window64 MS MPI V10.1.2, and install both msmpisdk.msi and MSMpiSetup.exe.  

In most cases just answer yes (or type return) to all questions.

Otherwise in the Additional Task windows, check the box "Add application directory to your system path." This is required otherwise the program ``ffglut.exe`` will not be found.

By now you should have two new icons on your desktop:

-  ``FreeFem++ (VERSION).exe``, the ``freefem++`` application.
-  ``FreeFem++ (VERSION) Examples``, a link to the ``freefem++`` examples folder.

where ``(VERSION)`` is the version of the files (for example 4.5).

By default, the installed files are in ``C:\Programs Files\FreeFem++``. In this directory, you have all the ``.dll`` files and other applications: ``FreeFem++-nw.exe``, ``ffglut.exe``, … The syntax for the command-line tools are the same as those of ``FreeFem.exe``.

To use FreeFEM binaries under Windows, two methods are possible:

* Use the FreeFEM launcher (launchff++.exe)

Warning: if you launch FreeFEM without filename script by double-clicking, your get a error due (it is bug of usage GetOpenFileName in win64).

* In shell terminal (cmd, powershell, bash, ... ):

  - To launch sequential version:   

  .. code-block:: bash
     :linenos:

     C:\>"Program Files (x86)\FreeFem++\FreeFem++.exe" <mySequentialScript.edp>
  
  - To launch parallel version:  
 
  .. code-block:: bash
     :linenos:

     C:\>"Program Files\Microsoft MPI\Bin\mpiexec.exe" -n <nbProcs> C:\>"Program Files (x86)\FreeFem++\FreeFem++-mpi.exe" <myParallelScript.edp>



macOS X installation
^^^^^^^^^^^^^^^^^^^^

Download the macOS X binary version file, extract all the files by double clicking on the icon of the file, go the the directory and put the ``FreeFem++.app`` application in the ``/Applications`` directory.

If you want terminal access to **FreeFEM** just copy the file ``FreeFem++`` in a directory of your :bash:`$PATH` shell environment variable.


Ubuntu installation
^^^^^^^^^^^^^^^^^^^

.. note:: The Debian package is built for Ubuntu 16.04   

Beforehand, install the following dependances libraries using the apt tool:

.. code-block:: bash
   :linenos:
    
   sudo apt-get install libgsl-dev libhdf5-dev 
                liblapack-dev libopenmpi-dev freeglut3-dev
	
Download the package FreeFEM .deb, install it by the command

.. code-block:: bash
   :linenos:
   
   dpkg -i FreeFEM_VERSION_Ubuntu_withPETSc_amd64.deb

FreeFEM is directly available in your terminal by the command "FreeFem++".


Arch AUR package
^^^^^^^^^^^^^^^^

An up-to-date package of **FreeFEM** for Arch is available on the `Archlinux user repository <https://aur.archlinux.org/packages/freefem%2B%2B-git/>`__.

To install it:

.. code-block:: bash
   :linenos:

   git clone https://aur.archlinux.org/freefem++-git.git
   cd freefem++-git
   makepkg -si

.. note:: Thanks to `Stephan Husmann <https://github.com/stefanhusmann>`__



Compiling source code
---------------------
	
Various versions of FreeFEM are possible: 
  - sequential and without plugins (contains in 3rdparty) 
  - parallel with plugins (and with PETSc).
  
   .. note:: We advise you to use the package manager for macOS Homebrew to get the different packages required avalaible `here <https://brew.sh>`__

Compilation on OSX (>=10.13)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Install Xcode, Xcode Command Line tools and Xcode Additional Tools from the `Apple website <https://developer.apple.com/download/more/>`__

2. Install gfortran from Homebrew

   .. code-block:: bash
      :linenos:

       brew install gfortran

3. To use **FreeFEM** parallel version, install `openmpi <https://www.open-mpi.org/software/ompi/v4.0/>`__  or  `mpich <http://www.mpich.org/downloads/>`__ 

   .. code-block:: bash
      :linenos:
	  
       # to install openmpi
       curl -L https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.1.tar.gz --output openmpi-4.0.1.tar.gz
       tar xf openmpi-4.0.1
       cd openmpi-4.0.1/
       # to install mpich
       curl -L http://www.mpich.org/static/downloads/3.3.2/mpich-3.3.2.tar.gz --output mpich-3.3.2.tar.gz
       tar xf mpich-3.3.2
       cd mpich-3.3.2
	   
   .. code-block:: bash
      :lineno-start: 4
	  
       # with brew gcc gfortran compilers
       ./configure CC=clang CXX=clang++ FC=gfortran-9 F77=gfortran-9 --prefix=/where/you/want/to/have/files/installed
    
       # with LLVM gcc and brew gfortran compilers
       ./configure CC=gcc-9 CXX=g++-9 FC=gfortran-9 F77=gfortran-9 --prefix=/where/you/want/to/have/files/installed
   

   .. code-block:: bash
      :lineno-start: 5

       make -j<nbProcs>
       make -j<nbProcs> install

4. Install the minimal libraries for **FreeFEM**

   .. code-block:: bash
      :linenos:

      brew install m4 git flex bison

5. If you want build your own configure according your system, install autoconf and automake from Homebrew (optional, see note in step 10)

   .. code-block:: bash
      :linenos:

      brew install autoconf automake

6. To use **FreeFEM** with its plugins, install from Homebrew suitesparse, hdf5, cmake, wget

   .. code-block:: bash
      :linenos:

      brew install suitesparse hdf5 cmake wget

7. Install `gsl <https://www.gnu.org/software/gsl>`__

   .. code-block:: bash
      :linenos:

      curl -O http://mirror.cyberbits.eu/gnu/gsl/gsl-2.5.tar.gz
      tar zxvf gsl-2.5.tar.gz
      cd gsl-2.5
      ./configure
      make -j<nbProcs>
      make -j<nbProcs> install --prefix=/where/you/want/to/have/files/installed

8. Download the latest Git for Mac installer `git <https://git-scm.com/download/mac>`__ and the **FreeFEM** source from the repository

   .. code-block:: bash
      :linenos:

      git clone https://github.com/FreeFem/FreeFem-sources.git

9. Configure your source code

   .. code-block:: bash
      :linenos:

       cd FreeFem-sources
       autoreconf -i
   .. note:: if your autoreconf version is too old, do ``tar zxvf AutoGeneratedFile.tar.gz``

   -  following your compilers

   .. code-block:: bash
      :lineno-start: 3
      
      // with brew gcc gfortran compilers 
      ./configure --enable-download CC=clang CXX=clang++ F77=gfortran-9 
	  FC=gfortran-9 --prefix=/where/you/want/to/have/files/installed
	  
      // with LLVM gcc and brew gfortran compilers 
      ./configure --enable-download CC=clang CXX=clang++ F77=gfortran-9 
	  FC=gfortran-9 --prefix=/where/you/want/to/have/files/installed

   

10. Download the 3rd party packages to use FreeFEM plugins

   .. code-block:: bash
      :linenos:

      ./3rdparty/getall -a

   .. note:: All the third party packages have their own licence

11. If you want use PETSc/SLEPc and `HPDDM <https://github.com/hpddm/hpddm>`__ (High Performance Domain Decomposition Methods)

   .. code-block:: bash
      :linenos:

      cd 3rdparty/ff-petsc
      make petsc-slepc // add SUDO=sudo if your installation directory is the default /usr/local
      cd -
      ./reconfigure

12. Build your **FreeFEM** library and executable

   .. code-block:: bash
      :linenos:

      make -j<nbProcs>
      make -j<nbProcs> check

   .. note:: ``make check`` is optionnally, but advise to check the validity of your **FreeFEM** building
   
13. Install the **FreeFEM** apllication 
      make install
     
   .. note:: it isn't necessary to execute this last command, FreeFEM executable is avalaible here your_installation/src/nw/FreeFem++ and mpi executable here your_installation/src/mpi/ff-mpirun



Compilation on Ubuntu
^^^^^^^^^^^^^^^^^^^^^

1. Install the following packages on your system

   .. code-block:: bash
      :linenos:

      sudo apt-get update && sudo apt-get upgrade
      sudo apt-get install cpp freeglut3-dev g++ gcc gfortran \
          m4 make patch pkg-config wget python unzip \
          liblapack-dev libhdf5-dev libgsl-dev \
          autoconf automake autotools-dev bison flex gdb git cmake

      # mpich is required for the FreeFEM parallel computing version
      sudo apt-get install mpich

   .. warning:: In the oldest distribution of Ubuntu, ``libgsl-dev`` does not exists, use ``libgsl2-dev`` instead

2. Download **FreeFEM** source from the repository

   .. code-block:: bash
      :linenos:

      git clone https://github.com/FreeFem/FreeFem-sources.git

3. Autoconf

   .. code-block:: bash
      :linenos:

      cd FreeFem-sources
      autoreconf -i

   .. note:: if your autoreconf version is too old, do ``tar zxvf AutoGeneratedFile.tar.gz``

4. Configure

   .. code-block:: bash
      :linenos:

      ./configure --enable-download --enable-optim 
	  --prefix=/where/you/want/to/have/files/installed

   .. note:: To see all the options, type ``./configure --help``

5. Download the 3rd party packages

   .. code-block:: bash
      :linenos:

      ./3rdparty/getall -a

   .. note:: All the third party packages have their own licence

6. If you want use PETSc/SLEPc and `HPDDM <https://github.com/hpddm/hpddm>`__ (High Performance Domain Decomposition Methods) for massively parallel computing

   .. code-block:: bash
      :linenos:

      cd 3rdparty/ff-petsc
      make petsc-slepc // add SUDO=sudo if your installation directory is the default /usr/local
      cd -
      ./reconfigure

7. Build your **FreeFEM** library and executable

   .. code-block:: bash
      :linenos:

      make -j<nbProcs>
      make -j<nbProcs> check
   
   .. note:: ``make check`` is optionnally, but advise to check the validity of your **FreeFEM** building

8. Install the executable 

   .. code-block:: bash
      :linenos:

      make install

   .. note:: it isn't necessary to execute this last command, FreeFEM executable is avalaible here your_installation/src/nw/FreeFem++ and mpi executable here your_installation/src/mpi/ff-mpirun


Compilation on Arch Linux
^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning:: As Arch is in rolling release, the following information can be quickly outdated !

.. warning:: **FreeFEM** fails to compile using the newest version of gcc 8.1.0, use an older one instead.

1. Install the following dependencies:

   .. code-block:: bash
      :linenos:

      pacman -Syu
      pacman -S git openmpi gcc-fortran wget python
          freeglut m4 make patch gmm
          blas lapack hdf5 gsl fftw arpack suitesparse
          gnuplot autoconf automake bison flex gdb
          valgrind cmake texlive-most

2. Download the **FreeFEM** source from the repository

   .. code-block:: bash
      :linenos:

      git clone https://github.com/FreeFem/FreeFem-sources.git

3. Autoconf

   .. code-block:: bash
      :linenos:

      cd FreeFem-sources
      autoreconf -i

4. Configure

   .. code-block:: bash
      :linenos:

      ./configure --enable-download --enable-optim

   .. note:: To see all the options, type ``./configure --help``

5. Download the packages

   .. code-block:: bash
      :linenos:

      ./3rdparty/getall -a

   .. note:: All the third party packages have their own licence

6. If you want use `HPDDM <https://github.com/hpddm/hpddm>`__ (High Performance Domain Decomposition Methods) for massively parallel computing, install PETSc/SLEPc

   .. code-block:: bash
      :linenos:

      cd 3rdparty/ff-petsc
      make petsc-slepc SUDO=sudo
      cd -
      ./reconfigure

7. Compile the **FreeFEM** source

   .. code-block:: bash
      :linenos:

      make

   .. note:: If your computer has many threads, you can run ``make`` in parallel using ``make -j16`` for 16 threads, for example.

   .. note:: Optionnally, check the compilation with ``make check``

8. Install the **FreeFEM** application

   .. code-block:: bash
      :linenos:

      sudo make install
	  
	  

Compilation on Linux with Intel software tools
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Follow the `guide <https://software.intel.com/en-us/articles/building-freefem-with-intel-software-tools-for-developers>`__


Compilation on Windows
^^^^^^^^^^^^^^^^^^^^^^

.. warning:: 
   The support ended for all releases under Windows 32 bits since the V4.
   We assume your development machine is 64-bit, and you want your compiler to target 64-bit windows by default.


1. Install the `Microsoft MPI v10.1.2 (archived) <https://www.microsoft.com/en-us/download/details.aspx?id=100593>`__ (msmpisdk.msi and MSMpiSetup.exe)

2. Download `msys2-x86_64-latest.exe <http://repo.msys2.org/distrib/msys2-x86_64-latest.exe>`__ (x86_64 version) and run it. 

3. Install the version control system `Git <https://git-scm.com/download/win>`__ for Windows

4. In the MSYS2 shell, execute the following. 
Hint: if you right click the title bar, go to Options -> Keys and tick "Ctrl+Shift+letter shortcuts" you can use Ctrl+Shift+V to paste in the MSYS shell.

   .. code-block:: bash
      :linenos:
      
      pacman -Syuu

Close the MSYS2 shell once you're asked to. There are now 3 MSYS subsystems installed: MSYS2, MinGW32 and MinGW64. 
They can respectively be launched from C:\dev\msys64\msys2.exe, C:\dev\msys64\mingw32.exe and C:\dev\msys64\mingw64.exe
Reopen MSYS2 (doesn't matter which version, since we're merely installing packages). 
Repeatedly run the following command until it says there are no further updates. You might have to restart your shell again.

   .. code-block:: bash
      :linenos:
      
      pacman -Syuu
	  

5. Now that MSYS2 is fully up-to-date, install the following dependancies

   -  for 64 bit systems:

   .. code-block:: bash
      :linenos:

      pacman -S autoconf make automake-wrapper bison git \
        mingw-w64-x86_64-freeglut mingw-w64-x86_64-toolchain \
        mingw-w64-x86_64-openblas patch python perl pkg-config pkgfile \
        rebase tar time tzcode unzip which mingw-w64-x86_64-libmicroutils \
        --ignore mingw-w64-x86_64-gcc-ada --ignore mingw-w64-x86_64-gcc-objc \
        --ignore mingw-w64-x86_64-gdb mingw-w64-x86_64-cmake --noconfirm

   -  for 32 bit systems (**FreeFEM** lower than version 4):

   .. code-block:: bash
      :linenos:

      pacman -S autoconf automake-wrapper bash bash-completion \
        bison bsdcpio bsdtar bzip2 coreutils curl dash file filesystem \
        findutils flex gawk gcc gcc-fortran gcc-libs grep gzip inetutils \
        info less lndir make man-db git mingw-w64-i686-freeglut \
        mingw-w64-i686-toolchain mingw-w64-i686-gsl mingw-w64-i686-hdf5 \
        mingw-w64-i686-openblas mintty msys2-keyring msys2-launcher-git \
        msys2-runtime ncurses pacman pacman-mirrors pactoys-git patch pax-git \
        perl pkg-config pkgfile rebase sed tar tftp-hpa time tzcode unzip \
        util-linux which

6. Open a ``MingW64 terminal`` (or ``MingW32`` for old 32 bit **FreeFEM** version) and compile the **FreeFEM** source

   .. code-block:: bash
      :linenos:

      git clone https://github.com/FreeFem/FreeFem-sources
      cd FreeFem-sources
      autoreconf -i
      ./configure --enable-generic --enable-optim \
        --enable-download --enable-maintainer-mode \
        CXXFLAGS=-mtune=generic CFLAGS=-mtune=generic \
        FFLAGS=-mtune=generic --enable-download --disable-hips
		--prefix=/where/you/want/to/have/files/installed
		
		
7. If you want use `HPDDM <https://github.com/hpddm/hpddm>`__ (High Performance Domain Decomposition Methods) for massively parallel computing, install PETSc/SLEPc

   .. code-block:: bash
      :linenos:

      cd 3rdparty/ff-petsc
      make petsc-slepc SUDO=sudo
      cd -
      ./reconfigure	

8. Download the 3rd party packages and build your **FreeFEM** library and executable

   .. code-block:: bash
      :linenos:		
	  
      ./3rdparty/getall -a
      make
      make check
      make install

   .. note:: The **FreeFEM** executable (and some other like ``ffmedit``, …) are in ``C:\msys64\mingw64\bin`` (or ``C:\msys32\mingw32\bin``).



.. .. _cmake:

.. Using CMake (FreeFEM without plugins)
   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. Compilation on OSX (>=10.13)
   """"""""""""""""""""""""""""

.. 1. Install Xcode, Xcode Command Line tools and Xcode Additional Tools from the `Apple website <https://developer.apple.com/download/more/>`__

.. 2. Install gcc from `http://hpc.sourceforge.net <http://hpc.sourceforge.net/>`__

..    .. code-block:: bash
..       :linenos:

..       curl -O http://prdownloads.sourceforge.net/hpc/gfortran-8.1-bin.tar.gz?download
..       sudo tar zxvf gfortran-8.1-bin.tar.gz -C /

.. 3. Install cmake from `macport <https://www.macports.org>`__ or with `Homebrew <https://brew.sh>`__

..    .. code-block:: bash
..       :linenos:

..       sudo port install cmake

..    .. code-block:: bash
..       :linenos:

..       brew install cmake

.. 4. Install mactex from `ctan <http://mirrors.ctan.org/systems/mac/mactex/MacTeX.pkg>`__

.. 5. Install the `openmpi <https://www.open-mpi.org/software/ompi/v4.0/>`__ source code

..    .. code-block:: bash
..       :linenos:

..       ./configure CC=/usr/local/bin/gcc CXX=/usr/local/bin/g++ F77=/usr/local/bin/gfortran FC=/usr/local/bin/gfortran
..       make
..       sudo make install

.. 6. Install `gsl <https://www.gnu.org/software/gsl>`__

..    .. code-block:: bash
..       :linenos:

..       curl -O https://fr.mirror.babylon.network/gnu/gsl/gsl-2.4.tar.gz
..       tar zxvf gsl-2.4.tar.gz
..       cd gsl-2.4
..       ./configure CC=/usr/local/bin/gcc
..       make
..       sudo make install

.. 7. Install `git <https://git-scm.com/download/mac>`__

.. 8. Install SparseSuite and Arpack from `macport <https://www.macports.org>`__ or with `Homebrew <https://brew.sh>`__

..   .. code-block:: bash
..       :linenos:

..       sudo port install arpack SuiteSparse

..    .. code-block:: bash
..       :linenos:

..       brew install arpack suite-sparse


.. 9. Download the **FreeFEM** source from the repository

..    .. code-block:: bash
..       :linenos:

..       git clone https://github.com/FreeFem/FreeFem-sources.git

.. 10. Compile **FreeFEM**. Don’t forget to update the MacOS SDK version with your own in the command below:

..    .. code-block:: bash
..       :linenos:

..       cd FreeFem-sources
..       mkdir build
..       cd build
..       cmake ..
..       make
..       make test
..       sudo make install


.. Compilation on Ubuntu
   """""""""""""""""""""

.. 1. Install the following dependencies

..    .. code-block:: bash
..       :linenos:

..       sudo apt-get update && sudo apt-get upgrade
..       sudo apt-get install cpp freeglut3-dev g++ gcc gfortran \
..           ghostscript m4 make patch pkg-config wget python unzip \
..           libopenblas-dev liblapack-dev libhdf5-dev libgsl-dev \
..           libscotch-dev libfftw3-dev libarpack2-dev libsuitesparse-dev \
..           libmumps-seq-dev libnlopt-dev coinor-libipopt-dev libgmm++-dev libtet1.5-dev \
..           gnuplot-qt autoconf automake autotools-dev bison flex gdb valgrind git cmake

..       # mpich is required for the FreeFem parallel computing version
..       sudo apt-get install mpich

..    .. warning:: In the oldest distribution of Ubuntu, ``libgsl-dev`` does not exists, use ``libgsl2-dev`` instead

.. 2. Download **FreeFEM** source from the repository

..    .. code-block:: bash
..       :linenos:

..       git clone https://github.com/FreeFem/FreeFem-sources.git

.. 3. Configure

..    .. code-block:: bash
..       :linenos:

..       cd FreeFem-sources
..       mkdir build
..       cd build
..       cmake ..

.. 4. Build

..   .. code-block:: bash
..       :linenos:

..       make

..    .. note:: If your computer has many threads, you can run ``make`` in parallel using ``make -j16`` for 16 threads, for example.

..    .. note:: Optionnally, check the compilation with ``make test``

.. 5. Install

..    .. code-block:: bash
..       :linenos:

..       sudo make install

..
.. Compilation on Windows
.. """"""""""""""""""""""
..
.. 1. Install `MS MPI v7 <https://www.microsoft.com/en-us/download/details.aspx?id=49926>`__ (msmpisdk.msi and MSMpiSetup.exe)
..
.. 2. Install `Msys2 <https://www.msys2.org/>`__ (x86_64 version)
..
.. 3. Start MSYS2 MSYS
..
.. 4. Open ``MSYS2 MSYS terminal`` to install dependancies
..
..    -  for 64bits system:
..
..    .. code-block:: bash
..       :linenos:
..
..       pacman -Syu
..       pacman -S autoconf automake-wrapper bash bash-completion \
..           bison bsdcpio bsdtar bzip2 cmake coreutils curl dash file filesystem \
..           findutils flex gawk gcc gcc-fortran gcc-libs grep gzip inetutils info less lndir \
..           make man-db git mingw-w64-x86_64-freeglut mingw-w64-x86_64-gcc \
..           mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-gsl mingw-w64-x86_64-hdf5 \
..           mingw-w64-x86_64-openblas mintty msys2-keyring msys2-launcher-git \
..           msys2-runtime ncurses pacman pacman-mirrors pactoys-git patch pax-git \
..           perl pkg-config pkgfile rebase sed tar tftp-hpa time tzcode unzip util-linux which
..
..    -  for 32bits system:
..
..    .. code-block:: bash
..       :linenos:
..
..       pacman -Syu
..       pacman -S autoconf automake-wrapper bash bash-completion \
..           bison bsdcpio bsdtar bzip2 cmake coreutils curl dash file filesystem \
..           findutils flex gawk gcc gcc-fortran gcc-libs grep gzip inetutils info less lndir \
..           make man-db git mingw-w64-i686-freeglut mingw-w64-i686-gcc \
..           mingw-w64-i686-gcc-fortran mingw-w64-i686-gsl mingw-w64-i686-hdf5 \
..           mingw-w64-i686-openblas mintty msys2-keyring msys2-launcher-git \
..           msys2-runtime ncurses pacman pacman-mirrors pactoys-git patch pax-git \
..       perl pkg-config pkgfile rebase sed tar tftp-hpa time tzcode unzip util-linux which
..
.. 5. Open ``MingW64 terminal`` (or ``MingW32``) to compile **FreeFEM**
..
..    .. code-block:: bash
..       :linenos:
..
..       git clone https://github.com/FreeFem/FreeFem-sources
..       cd FreeFem-sources
..       mkdir build
..       cd build
..       cmake ..
..       make -j4
..       make test
..       make install
..
..    The **FreeFEM** executable (and some other like ``ffmedit``, …)
..    are in ``C:\msys64\mingw64\bin`` (or ``C:\msys32\mingw32\bin``).





.. only:: html

  **FreeFEM** continuous integration 
  ----------------------------------

  The Inria Jenkins platform is used for the CI/CD integration of the source code.

  Compilation results of the develop branch by OS type and configuration of FreeFEM are here

  +------------------------+-------------------+-------------------+-------------------+-------------------+---------------------+
  | Branch                 | Linux 16.04       | Linux 18.04       | MacOS 10.10.5     | MacOS 10.13.5     | Windows 7           |
  +========================+===================+===================+===================+===================+=====================+
  | Develop                | |Build Status01|  | |Build Status02|  | |Build Status03|  | |Build Status04|  | |Build Status05|    |      
  +------------------------+-------------------+-------------------+-------------------+-------------------+---------------------+
  | Develop                | |Build Status06|  | |Build Status07|  | |Build Status08|  | |Build Status09|  | |Build Status10|    |
  | mpich and PETSc/SLEPc  |                   |                   |                   |                   | MSMPI V10.1.2       |
  +------------------------+-------------------+-------------------+-------------------+-------------------+---------------------+           
  | Develop                | |Build Status11|  | |Build Status12|  | |Build Status13|  | |Build Status14|  |                     |
  | openmpi and PETSc/SLEPc|                   |                   |                   |                   |                     |
  +------------------------+-------------------+-------------------+-------------------+-------------------+---------------------+ 


Environment variables and init file
-----------------------------------

**FreeFEM** reads a user’s init file named ``freefem++.pref`` to initialize global variables: :freefem:`verbosity`, :freefem:`includepath`, :freefem:`loadpath`.

.. note:: The variable :freefem:`verbosity` changes the level of internal printing (0: nothing unless there are syntax errors, 1: few, 10: lots, etc. …), the default value is 2.

   The included files are found in the :freefem:`includepath` list and the load files are found in the :freefem:`loadpath` list.

The syntax of the file is:

.. code-block:: bash
   :linenos:

   verbosity = 5
   loadpath += "/Library/FreeFem++/lib"
   loadpath += "/Users/hecht/Library/FreeFem++/lib"
   includepath += "/Library/FreeFem++/edp"
   includepath += "/Users/hecht/Library/FreeFem++/edp"
   # This is a comment
   load += "funcTemplate"
   load += "myfunction"
   load += "MUMPS_seq"

The possible paths for this file are

-  under Unix and MacOs

.. code-block:: bash
   :linenos:

   /etc/freefem++.pref
   $(HOME)/.freefem++.pref
   freefem++.pref

-  under windows

.. code-block:: bash
   :linenos:

   freefem++.pref

We can also use shell environment variables to change verbosity and the search rule before the init files.

.. code-block:: bash
   :linenos:

   export FF_VERBOSITY=50
   export FF_INCLUDEPATH="dir;;dir2"
   export FF_LOADPATH="dir;;dir3"

.. note:: The separator between directories must be ";" and not ":" because ":" is used under Windows.

.. note:: To show the list of init of **FreeFEM** , do

   .. code-block:: bash
      :linenos:

      export FF_VERBOSITY=100;
      ./FreeFem++-nw
	  
|
|  
	  

Coloring Syntax FreeFem++
-------------------------

Atom
^^^^

In order to get the syntax highlighting in `Atom <https://atom.io/>`__, you have to install the `FreeFEM language support <https://atom.io/packages/language-freefem-official>`__.

You can do it directly in Atom: Edit -> Preferences -> Install, and search for ``language-freefem-offical``.

To launch scripts directly from Atom, you have to install the ``atom-runner`` package. Once installed, modify the Atom configuration file (Edit -> Config...) to have something like that:

.. code-block:: bash
   :linenos:

   "*":
      ...

      runner:
         extensions:
            edp: "FreeFem++"
         scopes:
            "Freefem++": "FreeFem++"

Reboot Atom, and use Alt+R to run a FreeFem++ script.

Gedit
^^^^^

In order to get the syntax highlighting in Gedit, you have to downlaod the `Gedit parser <https://github.com/FreeFem/FreeFem-parser-gedit>`__ and copy it in ``/usr/share/gtksourceview-3.0/language-specs/``.	  


Textmate 2, an editor under macOS 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To use the coloring **FreeFEM** syntax with the Textmate 2 editor on Mac 10.7 or better, download from macromates.com and download the textmate freefem++ syntax `here <http://www3.freefem.org/ff++/Textmate2-ff++.zip>`__ (version june 2107). To install this parser, unzip Textmate2-ff++.zip and follow the explanation given in file How_To.rtf.

rom www.freefem.org/ff++/Textmate2-ff++.zip (version june 2107) unzip Textmate2-


Notepad++,an editor under windows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 
Read and follow the instruction, `FREEFEM++ COLOR SYNTAX OF WINDOWS <http://www3.freefem.org/ff++/color-syntax-win.pdf>`__ .
 
Emacs editor
^^^^^^^^^^^^
For emacs editor you can download `ff++-mode.el <https://github.com/rrgalvan/freefem-mode/>`__ .
 	  
	  
.. |Build Status01| image:: https://ci.inria.fr/freefem-dev/buildStatus/icon?job=FreeFEM-sources-ubuntu1604-job3
   :target: https://ci.inria.fr/freefem-dev/job/FreeFEM-sources-ubuntu1604-job3/
.. |Build Status02| image:: https://ci.inria.fr/freefem-dev/buildStatus/icon?job=FreeFEM-sources-ubuntu1804-job3
   :target: https://ci.inria.fr/freefem-dev/job/FreeFEM-sources-ubuntu1804-job3/
.. |Build Status03| image:: https://ci.inria.fr/freefem-dev/buildStatus/icon?job=FreeFEM-sources-macos1010-job3
   :target: https://ci.inria.fr/freefem-dev/job/FreeFEM-sources-macos1010-job3/
.. |Build Status04| image:: https://ci.inria.fr/freefem-dev/buildStatus/icon?job=FreeFEM-sources-macos1013-job3
   :target: https://ci.inria.fr/freefem-dev/job/FreeFEM-sources-macos1013-job3/
.. |Build Status05| image:: https://ci.inria.fr/freefem-dev/buildStatus/icon?job=FreeFEM-sources-windows7-job3
   :target: https://ci.inria.fr/freefem-dev/job/FreeFEM-sources-windows7-job3/

.. |Build Status06| image:: https://ci.inria.fr/freefem-dev/buildStatus/icon?job=FreeFEM-sources-ubuntu1604-job5_mpich/
   :target: https://ci.inria.fr/freefem-dev/job/FreeFEM-sources-ubuntu1604-job5/_mpich/
.. |Build Status07| image:: https://ci.inria.fr/freefem-dev/buildStatus/icon?job=FreeFEM-sources-ubuntu1804-job5_mpich/
   :target: https://ci.inria.fr/freefem-dev/job/FreeFEM-sources-ubuntu1804-job5/_mpich/
.. |Build Status08| image:: https://ci.inria.fr/freefem-dev/buildStatus/icon?job=FreeFEM-sources-macos1010-job5_mpich/
   :target: https://ci.inria.fr/freefem-dev/job/FreeFEM-sources-macos1010-job5/_mpich/
.. |Build Status09| image:: https://ci.inria.fr/freefem-dev/buildStatus/icon?job=FreeFEM-sources-macos1013-job5_mpich/
   :target: https://ci.inria.fr/freefem-dev/job/FreeFEM-sources-macos1013-job5/_mpich/
.. |Build Status10| image:: https://ci.inria.fr/freefem-dev/buildStatus/icon?job=FreeFEM-sources-windows7-job5
   :target: https://ci.inria.fr/freefem-dev/job/FreeFEM-sources-windows7-job5
   
.. |Build Status11| image:: https://ci.inria.fr/freefem-dev/buildStatus/icon?job=FreeFEM-sources-ubuntu1604-job5_openmpi/
   :target: https://ci.inria.fr/freefem-dev/job/FreeFEM-sources-ubuntu1604-job5/_openmpi/
.. |Build Status12| image:: https://ci.inria.fr/freefem-dev/buildStatus/icon?job=FreeFEM-sources-ubuntu1804-job5_openmpi/
   :target: https://ci.inria.fr/freefem-dev/job/FreeFEM-sources-ubuntu1804-job5/_openmpi/
.. |Build Status13| image:: https://ci.inria.fr/freefem-dev/buildStatus/icon?job=FreeFEM-sources-macos1010-job5_openmpi/
   :target: https://ci.inria.fr/freefem-dev/job/FreeFEM-sources-macos1010-job5/_openmpi/
.. |Build Status14| image:: https://ci.inria.fr/freefem-dev/buildStatus/icon?job=FreeFEM-sources-macos1013-job5_openmpi/
   :target: https://ci.inria.fr/freefem-dev/job/FreeFEM-sources-macos1013-job5/_openmpi/


