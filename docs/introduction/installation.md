# Installation guide

## Compilation on OSX (>=10.13)

Remark: All line beginning by # are shell command in terminal. 

1) Install Xcode, Xcode command line tools and Xcode Auxiliaire Tools from [Apple website](https://developer.apple.com/download/more/)

2) Install  gcc from [http://hpc.sourceforge.net](http://hpc.sourceforge.net/)

```
# curl -O http://prdownloads.sourceforge.net/hpc/gfortran-7.1-bin.tar.gz?download
# sudo tar zxvf gfortran-7.1-bin.tar.gz -C /
```

3) Install autoconf and automake from [macport](http://www.macports.org) or with [Homebrew](https://brew.sh)

```
# sudo port install  autoconf
# sudo port install  automatke
```

4) install mactex  from  [ctan](http://mirrors.ctan.org/systems/mac/mactex/MacTeX.pkg) 

5) install [openmpi](https://www.open-mpi.org/software/) source code

```
#  ./configure CC=/usr/local/bin/gcc CXX=/usr/local/bin/g++ F77=/usr/local/bin/gfortran FC=/usr/local/bin/gfortran
# make 
# sudo make install
```

6) Install [gsl](https://www.gnu.org/software/gsl) 

```
# curl -O https://fr.mirror.babylon.network/gnu/gsl/gsl-2.4.tar.gz
# tar zxvf gsl-2.4.tar.gz
# cd gsl-2.4
#./configure CC=/usr/local/bin/gcc
# make
# sudo make install 
```

7) Install [git](https://git-scm.com/download/mac)
    
8) Download FreeFem++ source from the repository

```
# git clone https://github.com/FreeFem/FreeFem-sources.git
```

9) Compile FreeFem++. Don't Forget to update the MacOSX sdk version with your own in the command below:

```
# cd ff++
# ./configure '-with-suffix=macos-10.13' '-without-fltk' '--enable-download' '--enable-optim' 'MPIRUN=/usr/local/bin/mpirun' '--enable-m64' '--without-x' 'CC=clang -isysroot /Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.13.sdk' 'CFLAGS=-mmacosx-version-min=10.13' 'CXXFLAGS=-mmacosx-version-min=10.13 -std=c++11' 'CXX=clang++ -isysroot /Applications/Xcode.app//Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.13.sdk' 'F77=/usr/local/bin/gfortran' 'FC=/usr/local/bin/gfortran' 'MPICXX=/usr/local/bin/mpic++' 'MPICC=/usr/local/bin/mpicc' 'MPIFC=/usr/local/bin/mpif90' 'MPIF77=/usr/local/bin/mpif90' '--enable-maintainer-mode'
# make
# sudo make install
```

## Compilation on Linux
TODO

## Compilation on Windows
TODO
