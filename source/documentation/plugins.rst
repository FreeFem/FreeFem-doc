.. role:: freefem(code)
  :language: freefem

.. role:: cpp(code)
 :language: cpp

Plugins
=======

gsl
---

The interface with ``gsl`` spline is available in **FreeFEM**, the seven kind of spline are

0. :freefem:`gslinterpcspline`: default type of spline
1. :freefem:`gslinterpakima`
2. :freefem:`gslinterpsteffen`
3. :freefem:`gslinterplinear`
4. :freefem:`gslinterppolynomial`
5. :freefem:`gslinterpcsplineperiodic`
6. :freefem:`gslinterpakimaperiodic`

A brief wing example given all the syntax:

.. code-block:: freefem
   :linenos:

   load "gsl"

   // Parameters
   int n = 10;
   real[int, int] dspline(2,n+1); //data points to define the spline
   for(int i = 0; i <= n; ++i){ //set data points
      real xx = square(real(i)/n);
      real yy = sin(xx*pi*2);
      dspline(0, i) = xx;
      dspline(1, i) = yy;
   }

   // GSL splines
   gslspline spline1(gslinterpcspline, dspline); //define the spline1
   gslspline spline11(dspline); //define the spline11
   gslspline spline2(gslinterpsteffen, dspline); //define the spline2
   gslspline spline3(gslinterpcspline, dspline(0, :), dspline(1, :));
   gslspline spline33(dspline(0, :), dspline(1, :)); //define the spline3
   spline1 = spline2; //copy spline2 in spline1

   real t = 1.;
   real s1 = spline1(t); //evaluate the function spline1 at t
   cout << "spline1(t) = " << s1 << endl;
   real ds1 = spline1.d(t); //evaluate the derivative of function spline1 at t
   cout << "spline1.d(t) = " << ds1 << endl;
   real dds1 = spline1.dd(t); //evaluate the second derivative of function spline1 at t
   cout << "spline1.dd(t) = " << dds1 << endl;

This can be usefull to build function from data value.

The list of all ``gsl`` functions and the **FreeFEM** equivalent is available in the :ref:`Language references <referenceFFGSLAWK>` (same names without ``_``).

ffrandom
--------

Plugin to linux ``random`` functions.

The range of the random generator is from :math:`0` to :math:`(2^{31})-1`.

.. code-block:: freefem
   :linenos:

   load "ffrandom"

   srandomdev(); //set a true random seed
   //warning: under window this command
   //change the seed by randinit(random())) so all
   //FreeFEM random function are changed

   int maxrang = 2^31 - 1;
   cout << " max range " << maxrang << endl;

   cout << random() << endl;
   cout << random() << endl;
   cout << random() << endl;

   srandom(10);
   cout << random() << endl;
   cout << random() << endl;
   cout << random() << endl;

mmap / semaphore
----------------

The idea is just try to use Interprocess communication using POSIX Shared Memory in Linux.

We build a small library ``libff-mmap-semaphore.c`` and ``libff-mmap-semaphore.h`` to easily interface.

-  mmap - allocate memory, or map files or devices into memory
-  semaphore - allow processes and threads to synchronize their actions

   A semaphore is an integer whose value is never allowed to fall below zero.
   Two operations can be performed on semaphores: increment the semaphore value by one (``sem_post``); and decrement the semaphore value by one (``sem_wait``).

   If the value of a semaphore is currently zero, then a ``sem_wait`` operation will block until the value becomes greater than zero.

**The functions of library**

First the ``semaphore`` interface to make synchronization:

-  :cpp:`typedef struct FF_P_sem *ff_Psem;` the pointer to data structure
-  :cpp:`ff_Psem ffsem_malloc();` malloc an empty data structure
-  :cpp:`void ffsem_del(ff_Psem sem);` clean and free the pointer
-  :cpp:`void ffsem_destroy(ff_Psem sem);` clean, close the data structure
-  :cpp:`void ffsem_init0(ff_Psem sem);` make a correct empty of the data structure
-  :cpp:`void ffsem_init(ff_Psem sem,const char *nmm, int crea);` create or use a new semaphore
-  :cpp:`long ffsem_post(ff_Psem sem);` ``nlocked``, the value of the semaphore is incremented, and all threads which are waiting on the semaphore are awakened
-  :cpp:`long ffsem_wait(ff_Psem sem);` the semaphore referenced by ``sem`` is locked.
   When calling ``sem_wait()``, if the semaphore’s value is zero, the calling thread will block until the lock is acquired or until the call is interrupted by a signal.

   Alternatively, the ``sem_trywait()`` function will fail if the semaphore is already locked, rather than blocking on the semaphore
-  :cpp:`long ffsem_trywait(ff_Psem p);`

Secondly, the ``mmap`` functions:

-  :cpp:`typedef struct FF_P_mmap *ff_Pmmap;` the pointer to data structure
-  :cpp:`ff_Psem ffmmap_malloc();` malloc an empty data structure
-  :cpp:`void ffmmap_del(ff_Pmmap p);` clean and free the pointer
-  :cpp:`void ffmmap_destroy(ff_Pmmap p);` clean, close the data structure
-  :cpp:`void ffmmap_init0(ff_Pmmap p);` make a correct empty of the data structure
-  :cpp:`long ffmmap_msync(ff_Pmmap p, long off, long ln);` call writes modified whole pages back to the filesystem and updates the file modification time.
   Only those pages containing ``addr`` and ``len-1`` succeeding locations will be examined.
-  :cpp:`void ffmmap_init(ff_Pmmap p, const char *nmm, long len);` allocate memory, or map files or devices into memory.
-  :cpp:`long ffmmap_read(ff_Pmmap p, void *t, size_t n, size_t off);` read ``n`` bytes from the ``mmap`` at memory ``off`` in pointer ``t``.
-  :cpp:`long ffmmap_write(ff_Pmmap p, void *t, size_t n, size_t off);` write ``n`` bytes to the ``mmap`` at memory ``off`` in pointer ``t``.

The **FreeFEM** corresponding functions:

-  :freefem:`Pmmap sharedata(filename, 1024);` new type to store the ``mmap`` informations of name store in string ``filename`` with 1024 is the size the ``sharedata`` zone and file.
-  :freefem:`Psemaphore smff("ff-slave", creat);` new type to store the semaphore of name ``ff-slave`` where ``creat`` is a boolean to create or use a existing semaphore.
-  :freefem:`Wait(sem)` the semaphore referenced by ``sem`` is locked.
   When calling :freefem:`Wait(sem)`, if the semaphore’s value is zero, the calling thread will block until the lock is acquired or until the call is interrupted by a signal.
   Alternatively, the :freefem:`trywait(sem)` function will fail if the semaphore is already locked, rather than blocking on the semaphore.
-  :freefem:`Post(sem)` the semaphore referenced by ``sem`` is unlocked, the value of the semaphore is incremented, and all threads which are waiting on the semaphore are awakened.
-  :freefem:`Read(sharedata ,offset, data);` read the variable ``data`` from the place ``offset`` in ``sharedata`` mmap.
-  :freefem:`Write(sharedata, offset, data);` write the variable ``data`` at the place ``offset`` in ``sharedata`` mmap.

The full example:

The ``FFMaster.c`` file:

.. code-block:: c
   :linenos:

   #include "libff-mmap-semaphore.h"
   #include <unistd.h>
   #include <stdlib.h>
   #include <stdio.h>
   ff_Psem sem_ff, sem_c; //the semaphore for mutex

   int main(int argc, const char ** argv)
   {
      int debug = 0;
      ff_Pmmap shd;
      double cff, rff;
      long status;
      int i;
      if (argc > 1) debug = atoi(argv[1]);
      ff_mmap_sem_verb = debug;

      sem_ff = ffsem_malloc();
      sem_c = ffsem_malloc();
      shd = ffmmap_malloc();

      ffsem_init(sem_ff, "ff-slave1", 1);
      ffsem_init(sem_c, "ff-master1", 1);
      ffmmap_init(shd, "shared-data", 1024);

      status = 1;
      ffmmap_write(shd, &status, sizeof(status), 8);
      ffmmap_msync(shd, 0, 32);

      char ff[1024];
      sprintf(ff, "FreeFem++ FFSlave.edp -nw -ns -v %d&", debug);
      system(ff); //lauch FF++ in batch no graphics
      if(debug) printf("cc: before wait\n");

      if(debug) printf("cc: before wait 0 ff\n");
      ffsem_wait(sem_ff);

      for (i = 0; i < 10; ++i){
         printf(" iter : %d \n", i);
         cff = 10+i;
         ffmmap_write(shd, &cff, sizeof(cff), 0);
         ffsem_post(sem_c);

         if(debug) printf(" cc: before wait 2\n");
         ffsem_wait(sem_ff);
         ffmmap_read(shd, &rff, sizeof(rff), 16);
         printf(" iter = %d rff= %f\n", i, rff);
      }

      status = 0; //end
      ffmmap_write(shd, &status, sizeof(status), 8);
      ffsem_post(sem_c);
      printf("End Master \n");
      ffsem_wait(sem_ff);
      ffsem_del(sem_ff);
      ffsem_del(sem_c);
      ffmmap_del(shd);
      return 0;
   }

The ``FFSlave.edp`` file:

.. code-block:: freefem
   :linenos:

   load "ff-mmap-semaphore"

   Psemaphore smff("ff-slave1", 0);
   Psemaphore smc("ff-master1", 0);
   Pmmap sharedata("shared-data", 1024);
   if (verbosity < 4) verbosity = 0;

   // Mesh
   mesh Th = square(10, 10);
   int[int] Lab = [1, 2, 3, 4];

   // Fespace
   fespace Vh(Th, P1);
   Vh u, v;

   // Macro
   macro grad(u) [dx(u), dy(u)] //

   int status = 1;
   cout << " FF status = " << status << endl;
   real cff, rff;

   // Problem
   problem Pb (u, v)
      = int2d(Th)(
           grad(u)'*grad(v)
      )
      - int2d(Th)(
           cff*v
      )
      + on(Lab, u=0)
      ;

   if (verbosity > 9) cout << " FF: before FF post\n";
   Post(smff); //unlock master end init

   while (1){
      if (verbosity > 9) cout << " FF: before FF wait \n";
      Wait(smc); //wait from cint write ok
      Read(sharedata, 0, cff);
      Read(sharedata, 8, status);

      cout << " After wait .. FF " << cff << " " << status << endl;
      if(status <= 0) break;

      // Solve
      Pb;
      rff = int2d(Th)(u*u);
      cout << " ** FF " << cff << " " << rff << endl;

      // Write
      Write(sharedata, 16, rff);
      Post(smff); //unlock cc
   }

   Post(smff); //wait from cint
   cout << " End FreeFEM " << endl;

To test this example of coupling ``C`` program and **FreeFEM** script:

.. code-block:: bash
   :linenos:

   cc -c libff-mmap-semaphore.c
   cc FFMaster.c -o FFMaster libff-mmap-semaphore.o -g -pthread
   ff-c++ -auto ff-mmap-semaphore.cpp
   ./FFMaster

The output:

.. code-block:: bash
   :linenos:

   len 1024 size 0
   len 1024 size 1024
   FF status = 1
   iter : 0
   After wait .. FF 10 1
   ** FF 10 0.161797
   iter = 0 rff= 0.161797
   iter : 1
   After wait .. FF 11 1
   ** FF 11 0.195774
   iter = 1 rff= 0.195774
   iter : 2
   After wait .. FF 12 1
   ** FF 12 0.232987
   iter = 2 rff= 0.232987
   iter : 3
   After wait .. FF 13 1
   ** FF 13 0.273436
   iter = 3 rff= 0.273436
   iter : 4
   After wait .. FF 14 1
   ** FF 14 0.317121
   iter = 4 rff= 0.317121
   iter : 5
   After wait .. FF 15 1
   ** FF 15 0.364042
   iter = 5 rff= 0.364042
   iter : 6
   After wait .. FF 16 1
   ** FF 16 0.414199
   iter = 6 rff= 0.414199
   iter : 7
   After wait .. FF 17 1
   ** FF 17 0.467592
   iter = 7 rff= 0.467592
   iter : 8
   After wait .. FF 18 1
   ** FF 18 0.524221
   iter = 8 rff= 0.524221
   iter : 9
   After wait .. FF 19 1
   ** FF 19 0.584086
   iter = 9 rff= 0.584086
   End Master
   After wait .. FF 19 0
