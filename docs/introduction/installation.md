# Installation guide

## Easy installation

First, open the following web page of [__`FreeFem++`__](http://www.freefem.org/).

Choose your platform: Linux, MacOS or Windows in the download section.

!!!note
	Binary packages are available for Microsoft Windows, Apple MacOS X and some Linux distributions.

Install by double-click on the appropriate file. Under Linux and MacOS the install directory is one of the following `/usr/local/bin`, `/usr/local/share/freefem++`, `/usr/local/lib/ff++`

### Windows binary installation

First download the windows installation executable, then double click on it to install __`FreeFem++`__.

In most cases just answer yes (or type return) to all questions.

Otherwise in the Additional Task windows, check the box ”Add application directory to your system path.” This is required otherwise the program `ffglut.exe` will not be found.

By now you should have two new icons on your desktop:

 * `FreeFem++ (VERSION).exe`, the __`FreeFem++`__ application.
 * `FreeFem++ (VERSION) Examples`, a link to the __`FreeFem++`__ folder of examples.

where `(VERSION)` is the version of the files (for example 3.59).

By default, the installed files are in `C:\Programs Files\FreeFem++`.

In this directory, you have all the `.dll` files and other applications: `FreeFem++-nw.exe`, `ffglut.exe`, ...

The syntax for the command-line tools are the same as those of `FreeFem.exe`.

### MacOS X binary installation

Download the MacOS X binary version file, extract all the files with a double click on the icon of the file, go the the directory and put the `FreeFem+.app` application in the `/Applications` directory.

If you want a terminal access to __`FreeFem++`__ just copy the file `FreeFem++` in a directory of your `$PATH` shell environment variable.

If you want to automatically launch the `FreeFem++.app`, double click on a `.edp` file icon.

Under the finder pick a `.edp`, select menu File -> Get Info and change Open with: (choose `FreeFem++.app`) and click
on button change All.

## Compilation

### Compilation on OSX (>=10.13)

Remark: Blocks of code are shell commands in terminal.

1) Install Xcode, Xcode command line tools and Xcode Auxiliaire Tools from [Apple website](https://developer.apple.com/download/more/)

2) Install gcc from [http://hpc.sourceforge.net](http://hpc.sourceforge.net/)

```
curl -O http://prdownloads.sourceforge.net/hpc/gfortran-7.1-bin.tar.gz?download
sudo tar zxvf gfortran-7.1-bin.tar.gz -C /
```

3) Install autoconf and automake from [macport](http://www.macports.org) or with [Homebrew](https://brew.sh)

```
sudo port install autoconf
sudo port install automake
```

4) install mactex from [ctan](http://mirrors.ctan.org/systems/mac/mactex/MacTeX.pkg)

5) install [openmpi](https://www.open-mpi.org/software/) source code

```
./configure CC=/usr/local/bin/gcc CXX=/usr/local/bin/g++ F77=/usr/local/bin/gfortran FC=/usr/local/bin/gfortran
make
sudo make install
```

6) Install [gsl](https://www.gnu.org/software/gsl)

```
curl -O https://fr.mirror.babylon.network/gnu/gsl/gsl-2.4.tar.gz
tar zxvf gsl-2.4.tar.gz
cd gsl-2.4
./configure CC=/usr/local/bin/gcc
make
sudo make install
```

7) Install [git](https://git-scm.com/download/mac)

8) Download __`FreeFem++`__ source from the repository

```
git clone https://github.com/FreeFem/FreeFem-sources.git
```

9) Compile __`FreeFem++`__. Don't Forget to update the MacOSX sdk version with your own in the command below:

```bash
cd FreeFem-sources
./configure '-with-suffix=macos-10.13' '-without-fltk' '--enable-download' '--enable-optim' 'MPIRUN=/usr/local/bin/mpirun' '--enable-m64' '--without-x' 'CC=clang -isysroot /Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.13.sdk' 'CFLAGS=-mmacosx-version-min=10.13' 'CXXFLAGS=-mmacosx-version-min=10.13 -std=c++11' 'CXX=clang++ -isysroot /Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.13.sdk' 'F77=/usr/local/bin/gfortran' 'FC=/usr/local/bin/gfortran' 'MPICXX=/usr/local/bin/mpic++' 'MPICC=/usr/local/bin/mpicc' 'MPIFC=/usr/local/bin/mpif90' 'MPIF77=/usr/local/bin/mpif90' '--enable-maintainer-mode'
make
sudo make install
```

### Compilation on Ubuntu

1) Install the following dependencies

```bash
sudo apt-get update && sudo apt-get upgrade
sudo apt-get install cpp freeglut3-dev g++ gcc gfortran \
	ghostscript m4 make patch pkg-config wget python unzip \
	libopenblas-dev liblapack-dev libhdf5-dev libgsl2-dev \
	libscotch-dev libfftw3-dev libarpack2-dev libsuitesparse-dev \
	libmumps-seq-dev libnlopt-dev coinor-libipopt-dev libgmm++-dev libtet1.5-dev \
	gnuplot-qt autoconf automake autotools-dev bison flex gdb valgrind git cmake

# mpich is required for the FreeFem parallel computing version
sudo apt-get install mpich
```

!!!warning
	In the latest distribution of Ubuntu, `libgsl2-dev` does not exists anymore, use `libgsl-dev`

2) Download __`FreeFem++`__ source from the repository

```bash
git clone https://github.com/FreeFem/FreeFem-sources.git
```

3) Autoconf

```bash
cd FreeFem-sources
autoreconf -i
```

!!!info
	if your autoreconf version is too old, do `tar zxvf AutoGeneratedFile.tar.gz`

4) Configure

```bash
./configure --enable-download --enable-optim --disable-pastix
```

!!!info
	To see all the options, type `./configure --help`

5) Download the packages

```bash
./download/getall -a
```

!!!info
	All the third party packages have their own licence

6) Download and compile petsc & slepc

```bash
cd download/ff-petsc
make petsc-slepc SUDO=sudo
cd -
```

!!! warning
	Pastix seems to fail during the PETSc compilation. You can remove `--download-pastix` directly in the `Makefile` if a compilation error occurs.

7) Reconfigure with petsc and slepc

```bash
./reconfigure
```

8) Build

```bash
make
```

!!!note
	If your computer have many threads, you can run `make` in parallel using `make -j16` for 16 threads for example.

!!!info
	Optionnally, check the compilation with `make check`

9) Install

```bash
sudo make install
```

### Compilation on Arch Linux

!!! warning
	As Arch is in rolling release, the following informations can be quickly exceeded !

1) Install the following dependencies:

```bash
pacman -Syu
pacman -S git openmpi gcc-fortran wget python
	freeglut ghostscript m4 make patch gmm
	blas lapack hdf5 gsl fftw arpack suitesparse
	gnuplot autoconf automake bison flex gdb
	valgrind cmake texlive-most

```

2) Download __`FreeFem++`__ source from the repository

```bash
git clone https://github.com/FreeFem/FreeFem-sources.git
```

3) Autoconf

```bash
cd FreeFem-sources
autoreconf -i
```

4) Configure

```bash
./configure --enable-download --enable-optim --disable-pastix
```

!!!info
	To see all the options, type `./configure --help`

5) Download the packages

```bash
./download/getall -a
```

!!!info
	All the third party packages have their own licence

6) Download and compile petsc & slepc

```bash
cd download/ff-petsc
make petsc-slepc SUDO=sudo
cd -
```
7) Reconfigure with petsc and slepc

```bash
./reconfigure
```

8) Build

```bash
make
```

!!!note
	If your computer have many threads, you can run `make` in parallel using `make -j16` for 16 threads for example.

!!!info
	Optionnally, check the compilation with `make check`

9) Install

```bash
sudo make install
```

### Compilation on Linux with Intel software tools

Follow the [guide](https://software.intel.com/en-us/articles/building-freefem-with-intel-software-tools-for-developers)

### Compilation on Windows
$\codered$ (Good luck!)

## Environnement variables and init file

__`FreeFem++`__ reads a user’s init file named `freefem++.pref` to initialize global variables: `:::freefem verbosity`, `:::freefem includepath`, `:::freefem loadpath`.

!!!note
	The variable `:::freefem verbosity` changes the level of internal printing (0, nothing unless there are syntax errors, 1 few, 10 lots, etc. ...), the default value is 2.

	The include files are searched from the `:::freefem includepath` list and the load files are searched from `:::freefem loadpath` list.

The syntax of the file is:

```bash
verbosity = 5
loadpath += "/Library/FreeFem++/lib"
loadpath += "/Users/hecht/Library/FreeFem++/lib"
includepath += "/Library/FreeFem++/edp"
includepath += "/Users/hecht/Library/FreeFem++/edp"
# a comment
load += "funcTemplate"
load += "myfunction"
load += "MUMPS_seq"
```

The possible paths for this file are

 * under Unix and MacOs
 ```bash
 /etc/freefem++.pref
 $(HOME)/.freefem++.pref
 freefem++.pref
 ```
* under windows
 ```bash
 freefem++.pref
 ```

We can also use shell environment variable to change verbosity and the search rule before the init files.
```bash
export FF_VERBOSITY=50
export FF_INCLUDEPATH="dir;;dir2"
export FF_LOADPATH="dir;;dir3"
```

!!!note
	The separator between directories must be ”;” and not ”:” because ”:” is used under Windows.

!!!note
	To show the list of init of FreeFem++ , do
	```bash
	export FF_VERBOSITY=100;
	./FreeFem++-nw
	```
