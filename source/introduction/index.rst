.. raw:: html

   <style type="text/css">
   t {
     color: #00007f;
     font-weight: bold;
	 font-size:50px
   }
   </style>

   <center>
    <t> FreeFEM </t> <br/>
    Fourth edition, Version 4.2.1 </t>
   </center>

.. only:: html

  .. raw:: html

    <div id="PDFbanner">
	  <p>
	    <strong><a href="https://www.ljll.math.upmc.fr/hecht/" target="_blank">Frédéric Hecht</a></strong>
    </p>
    <p>
		  Laboratoire Jacques-Louis Lions, Sorbonne University, Paris
	  </p>
    <p>
      The user manual is available <a href="/pdf/FreeFEM-documentation.pdf"><img src="../_static/img/logo_pdf.png"></a>
    </p>
	  <div class="PDFbanner-v3">
        <p>
          <i>The FreeFEM v3 PDF is archived (not up to date)</i>
        </p>
        <p>
          <a href="/_static/pdf/FreeFEM-doc-v3.pdf" target="_blank">Offical version</a> (English)<br/>
          <a href="/_static/pdf/FreeFEM-doc-v3_Spanish.pdf" target="_blank">Spanish version</a>, by Eliseo Chac ́on Vera <br/>
          <a href="/_static/pdf/FreeFEM-doc-v3_Chinese.pdf" target="_blank">Chinese version</a> (中文) <br/>
		      <a href="http://comfos.org/jp/ffempp/index.html" target="_blank">Japanese version</a>, by Kohji Ohtsuka
        </p>
      </div>
    </div>

**Introduction**
================

**FreeFEM** is a partial differential equation solver for non-linear multi-physics systems in 2D and 3D.

Problems involving partial differential equations from several branches of physics, such as fluid-structure interactions, require interpolations of data on several meshes and their manipulation within one program.
**FreeFEM** includes a fast interpolation algorithm and a language for the manipulation of data on multiple meshes.

**FreeFEM** is written in C++ and its language is a C++ idiom.

**FreeFEM** currently interfaces to the following libraries:

.. hlist::
  :columns: 3

  * `ARPACK <https://www.caam.rice.edu/software/ARPACK/>`_
  * `BLAS <http://www.netlib.org/blas/>`_
  * `OpenBLAS <http://www.openblas.net/>`_
  * `FFTW V3.3.2 <http://www.fftw.org>`_
  * `Ipopt V3.12.4 <https://github.com/coin-or/Ipopt>`_
  * `Gmm++ V4.2 <http://getfem.org/gmm.html>`_
  * `freeYams <https://www.ljll.math.upmc.fr/frey/software.html>`_
  * `METIS V5.1.0 <http://glaros.dtc.umn.edu/gkhome/metis/metis/overview>`_
  * `ParMETIS V4.0.3 <http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview>`_
  * `MMG3D V4.0 <https://www.mmgtools.org>`_
  * `mshmet <https://www.ljll.math.upmc.fr/frey/software.html>`_
  * `MUMPS V5.0.2 <http://mumps.enseeiht.fr/>`_
  * `NLopt V2.2.4 <http://ab-initio.mit.edu/wiki/index.php/NLopt>`_
  * `ScaLAPACK <http://www.netlib.org/scalapack/>`_
  * `Scotch V6.0.4 <https://gforge.inria.fr/projects/scotch/>`_
  * `SuiteSparse V4.4.4 <http://faculty.cse.tamu.edu/davis/suitesparse.html>`_
  * `SuperLU V5.2.1 <http://crd-legacy.lbl.gov/~xiaoye/SuperLU/>`_
  * `TetGen V1.5.1 <http://www.tetgen.org/>`_
  * `PETSc V3.11.2 <https://www.mcs.anl.gov/petsc/>`_
  * `HTool <https://pierremarchand.netlify.com/project/htool/>`_
  * `HPDDM <https://github.com/hpddm/hpddm>`_

.. _FreeFem++: https://freefem.org

.. only:: html

  .. raw:: html

    <div id="IndexLogos">
      <img src="../_static/img/logo_cnrs_SU_UP.png" alt="CNRS"/>
      <img src="../_static/img/inr_logo_rouge_300.jpg" alt="INRIA"/>
      <img src="../_static/img/logo_LJLL.png" alt="LJLL"/>
      <img src="../_static/img/logo_alpines.png" alt="Alpines"/>
    </div>

.. toctree::

   new-features
   installation
   download
   history
   citation
   authors
   contributing
