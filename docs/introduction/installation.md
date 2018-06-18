# Installation guide

## Easy installation

First, go to the [__download__](download.md) page and choose your platform: Linux, MacOS or Windows.

!!!note
	Binary packages are available for Microsoft Windows, MacOS and some Linux distributions.

Install `FreeFem++` by double-clicking on the appropriate file. Under Linux and MacOS the install directory is one of the following `/usr/local/bin`, `/usr/local/share/freefem++`, `/usr/local/lib/ff++`

### Windows binary installation

First download the windows installation executable, then double click to install `freefem++`.

In most cases just answer yes (or type return) to all questions.

Otherwise in the Additional Task windows, check the box ”Add application directory to your system path.” This is required otherwise the program `ffglut.exe` will not be found.

By now you should have two new icons on your desktop:

 * `FreeFem++ (VERSION).exe`, the `freefem++` application.
 * `FreeFem++ (VERSION) Examples`, a link to the `freefem++` examples folder.

where `(VERSION)` is the version of the files (for example 3.59).

By default, the installed files are in `C:\Programs Files\FreeFem++`. In this directory, you have all the `.dll` files and other applications: `FreeFem++-nw.exe`, `ffglut.exe`, ...
The syntax for the command-line tools are the same as those of `FreeFem.exe`.

### MacOS X binary installation

Download the MacOS X binary version file, extract all the files by double clicking on the icon of the file, go the the directory and put the `FreeFem+.app` application in the `/Applications` directory.

If you want terminal access to `freefem++` just copy the file `FreeFem++` in a directory of your `$PATH` shell environment variable.

## Compilation

### Branches / OS status

| Branch | Linux | MacOSX | Windows 7 |
|:---:|:---:|:---:|:---:|
| Develop | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-develop-UbuntuAll)](https://ci.inria.fr/freefem/job/FreeFem-source-develop-UbuntuAll/) [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-develop-UbuntuNo)](https://ci.inria.fr/freefem/job/FreeFem-source-develop-UbuntuNo/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-develop-MacOSXAll)](https://ci.inria.fr/freefem/job/FreeFem-source-develop-MacOSXAll/) [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-develop-MacOSXNo)](https://ci.inria.fr/freefem/job/FreeFem-source-develop-MacOSXNo/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-develop-Windows7)](https://ci.inria.fr/freefem/job/FreeFem-source-develop-Windows7) [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-develop-Windows7-32)](https://ci.inria.fr/freefem/job/FreeFem-source-develop-Windows7-32) |
| Master | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-master-UbuntuAll)](https://ci.inria.fr/freefem/job/FreeFem-source-master-UbuntuAll/) [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-master-UbuntuNo)](https://ci.inria.fr/freefem/job/FreeFem-source-master-UbuntuNo/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-master-MacOSXAll)](https://ci.inria.fr/freefem/job/FreeFem-source-master-MacOSXAll/) [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-master-MacOSXNo)](https://ci.inria.fr/freefem/job/FreeFem-source-master-MacOSXNo/) | [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-master-Windows7)](https://ci.inria.fr/freefem/job/FreeFem-source-master-Windows7) [![Build Status](https://ci.inria.fr/freefem/buildStatus/icon?job=FreeFem-source-master-Windows7-32)](https://ci.inria.fr/freefem/job/FreeFem-source-master-Windows7-32) |

### Compilation on OSX (>=10.13)

Remark: Blocks of code are shell commands in terminal.

1. Install Xcode, Xcode Command Line tools and Xcode Additional Tools from the [Apple website](https://developer.apple.com/download/more/)

2. Install gcc from [http://hpc.sourceforge.net](http://hpc.sourceforge.net/)

	```bash
	curl -O http://prdownloads.sourceforge.net/hpc/gfortran-7.1-bin.tar.gz?download
	sudo tar zxvf gfortran-7.1-bin.tar.gz -C /
	```

3. Install autoconf and automake from [macport](http://www.macports.org) or with [Homebrew](https://brew.sh)

	```bash
	sudo port install autoconf
	sudo port install automake
	```

4. Install mactex from [ctan](http://mirrors.ctan.org/systems/mac/mactex/MacTeX.pkg)

5. Install the [openmpi](https://www.open-mpi.org/software/) source code

	```bash
	./configure CC=/usr/local/bin/gcc CXX=/usr/local/bin/g++ F77=/usr/local/bin/gfortran FC=/usr/local/bin/gfortran
	make
	sudo make install
	```

6. Install [gsl](https://www.gnu.org/software/gsl)

	```bash
	curl -O https://fr.mirror.babylon.network/gnu/gsl/gsl-2.4.tar.gz
	tar zxvf gsl-2.4.tar.gz
	cd gsl-2.4
	./configure CC=/usr/local/bin/gcc
	make
	sudo make install
	```

7. Install [git](https://git-scm.com/download/mac)

8. Download the FreeFem++ source from the repository

	```bash
	git clone https://github.com/FreeFem/FreeFem-sources.git
	```

9) Compile FreeFem++. Don't forget to update the MacOS SDK version with your own in the command below:

	```bash
	cd FreeFem-sources
	./configure '-with-suffix=macos-10.13' '-without-fltk' '--enable-download' '--enable-optim' 'MPIRUN=/usr/local/bin/mpirun' '--enable-m64' '--without-x' 'CC=clang -isysroot /Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.13.sdk' 'CFLAGS=-mmacosx-version-min=10.13' 'CXXFLAGS=-mmacosx-version-min=10.13 -std=c++11' 'CXX=clang++ -isysroot /Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.13.sdk' 'F77=/usr/local/bin/gfortran' 'FC=/usr/local/bin/gfortran' 'MPICXX=/usr/local/bin/mpic++' 'MPICC=/usr/local/bin/mpicc' 'MPIFC=/usr/local/bin/mpif90' 'MPIF77=/usr/local/bin/mpif90' '--enable-maintainer-mode'
	make
	sudo make install
	```

### Compilation on Ubuntu

1. Install the following dependencies

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

2. Download `freefem++` source from the repository

	```bash
	git clone https://github.com/FreeFem/FreeFem-sources.git
	```

3. Autoconf

	```bash
	cd FreeFem-sources
	autoreconf -i
	```

	!!!info
		if your autoreconf version is too old, do `tar zxvf AutoGeneratedFile.tar.gz`

4. Configure

	```bash
	./configure --enable-download --enable-optim --disable-pastix
	```

	!!!info
		To see all the options, type `./configure --help`

5. Download the packages

	```bash
	./download/getall -a
	```

	!!!info
		All the third party packages have their own licence

6. Download and compile petsc & slepc

	```bash
	cd download/ff-petsc
	make petsc-slepc SUDO=sudo
	cd -
	```

	!!! warning
		Pastix seems to fail during the PETSc compilation. You can remove `--download-pastix` directly in the `Makefile` if a compilation error occurs.

7. Reconfigure with petsc and slepc

	```bash
	./reconfigure
	```

8. Build

	```bash
	make
	```

	!!!note
		If your computer has many threads, you can run `make` in parallel using `make -j16` for 16 threads, for example.

	!!!info
		Optionnally, check the compilation with `make check`

9. Install

	```bash
	sudo make install
	```

### Compilation on Arch Linux

!!! warning
	As Arch is in rolling release, the following information can be quickly outdated !

!!!warning
	**`FreeFem++`** fails to compile using the newest version of gcc 8.1.0, use an older one instead.

1. Install the following dependencies:

	```bash
	pacman -Syu
	pacman -S git openmpi gcc-fortran wget python
		freeglut ghostscript m4 make patch gmm
		blas lapack hdf5 gsl fftw arpack suitesparse
		gnuplot autoconf automake bison flex gdb
		valgrind cmake texlive-most
	```

2. Download the `FreeFem++` source from the repository

	```bash
	git clone https://github.com/FreeFem/FreeFem-sources.git
	```

3. Autoconf

	```bash
	cd FreeFem-sources
	autoreconf -i
	```

4. Configure

	```bash
	./configure --enable-download --enable-optim --disable-pastix
	```

	!!!info
		To see all the options, type `./configure --help`

5. Download the packages

	```bash
	./download/getall -a
	```

	!!!info
		All the third party packages have their own licence

6. Download and compile petsc & slepc

	```bash
	cd download/ff-petsc
	make petsc-slepc SUDO=sudo
	cd -
	```

7. Reconfigure with petsc and slepc

	```bash
	./reconfigure
	```

8. Build

	```bash
	make
	```

	!!!note
		If your computer has many threads, you can run `make` in parallel using `make -j16` for 16 threads, for example.

	!!!info
		Optionnally, check the compilation with `make check`

9. Install

	```bash
	sudo make install
	```

### Compilation on Linux with Intel software tools

Follow the [guide](https://software.intel.com/en-us/articles/building-freefem-with-intel-software-tools-for-developers)

### Compilation on Windows

1. Install [MS MPI v7](https://www.microsoft.com/en-us/download/details.aspx?id=49926) (msmpisdk.msi and MSMpiSetup.exe)

2. Install [Msys2](https://www.msys2.org/) (x86_64 version)

3. Start MSYS2 MSYS

4. Open `MSYS2 MSYS terminal` to install dependancies
	- for 64bits system:
	```bash
	pacman -Syu
	pacman -S autoconf automake-wrapper bash bash-completion \
	bison bsdcpio bsdtar bzip2 coreutils curl dash file filesystem \
	findutils flex gawk gcc gcc-fortran gcc-libs grep gzip inetutils info less lndir \
	make man-db git mingw-w64-x86_64-freeglut mingw-w64-x86_64-gcc \
	mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-gsl mingw-w64-x86_64-hdf5 \
	mingw-w64-x86_64-openblas mintty msys2-keyring msys2-launcher-git \
	msys2-runtime ncurses pacman pacman-mirrors pactoys-git patch pax-git \
	perl pkg-config pkgfile rebase sed tar tftp-hpa time tzcode unzip util-linux which
	```
	- for 32bits system:
	```bash
	pacman -Syu
	pacman -S autoconf automake-wrapper bash bash-completion \
	bison bsdcpio bsdtar bzip2 coreutils curl dash file filesystem \
	findutils flex gawk gcc gcc-fortran gcc-libs grep gzip inetutils info less lndir \
	make man-db git mingw-w64-i686-freeglut mingw-w64-i686-gcc \
	mingw-w64-i686-gcc-fortran mingw-w64-i686-gsl mingw-w64-i686-hdf5 \
	mingw-w64-i686-openblas mintty msys2-keyring msys2-launcher-git \
	msys2-runtime ncurses pacman pacman-mirrors pactoys-git patch pax-git \
	perl pkg-config pkgfile rebase sed tar tftp-hpa time tzcode unzip util-linux which
	```
5. Open `MingW64 terminal` (or `MingW32`) to compile __`FreeFem++`__
	```bash
	git clone https://github.com/FreeFem/FreeFem-sources
	cd FreeFem-sources
	autoreconf -i
	./configure --enable-download --disable-pastix --disable-hips
	./download/getall -a
	make -j4
	make check
	make install
	```

	The __`FreeFem++`__ executable (and some other like `ffmedit`, ...) are in `C:\msys64\mingw64\bin` (or `C:\msys32\mingw64\bin`).

## Environment variables and init file

__`FreeFem++`__ reads a user’s init file named `freefem++.pref` to initialize global variables: `:::freefem verbosity`, `:::freefem includepath`, `:::freefem loadpath`.

!!!note
	The variable `:::freefem verbosity` changes the level of internal printing (0: nothing unless there are syntax errors, 1: few, 10: lots, etc. ...), the default value is 2.

	The included files are found in the `:::freefem includepath` list and the load files are found in the `:::freefem loadpath` list.

The syntax of the file is:

```bash
verbosity = 5
loadpath += "/Library/FreeFem++/lib"
loadpath += "/Users/hecht/Library/FreeFem++/lib"
includepath += "/Library/FreeFem++/edp"
includepath += "/Users/hecht/Library/FreeFem++/edp"
# This is a comment
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

We can also use shell environment variables to change verbosity and the search rule before the init files.

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
