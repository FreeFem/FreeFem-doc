|
|


.. raw:: html
   
   <style type="text/css">
   t {
     color: blue;
     font-weight: bold;
	 font-size:50px
   }
   </style>
   
   <center>
    <t> FreeFEM </t> <br/>
    Fourth edition, Version 4.8 </t>
   </center> 

===============

|
	
.. only:: html

  .. raw:: html

    <div id="PDFbanner">
	  <p>
	    <b><font size="5"><a href="https://www.ljll.math.upmc.fr/hecht/" target="_blank">Frédéric Hecht</a></font></b>  <br> 
	    
		<font size="4">Laboratoire Jacques-Louis Lions, Sorbonne University, Paris</font>
	  </p>	
	  <br> 
      <p>
        <font size="4"> The user manual is available <a href="/pdf/FreeFEM-documentation.pdf"><img src="../_static/img/logo_pdf.png"width="25" height=""></a></font>
      </p>  
	  <br> 
	  <br> 
	  <div class="PDFbanner-v3">
        <p>
          <i>The FreeFEM v3 PDF is archived (not up to date)</i>
        </p>
        <p>
          <a href="/_static/pdf/FreeFEM-doc-v3.pdf" target="_blank">Official version</a> (English)<br/>
          <a href="/_static/pdf/FreeFEM-doc-v3_Spanish.pdf" target="_blank">Spanish version</a>, by Eliseo Chacón Vera <br/>
          <a href="/_static/pdf/FreeFEM-doc-v3_Chinese.pdf" target="_blank">Chinese version</a> (中文) <br/>
		  <a href="http://comfos.org/jp/ffempp/index.html" target="_blank">Japanese version</a>, by Kohji Ohtsuka
        </p>
      </div>
    </div>
	
|
|
|
|

**Introduction**
================

**FreeFEM** is a partial differential equation solver for non-linear multi-physics systems in 1D, 2D, 3D and 3D border domains (surface and curve).

Problems involving partial differential equations from several branches of physics, such as fluid-structure interactions, require interpolations of data on several meshes and their manipulation within one program.
**FreeFEM** includes a fast interpolation algorithm and a language for the manipulation of data on multiple meshes.

**FreeFEM** is written in C++ and its language is a C++ idiom.

|
|

**FreeFEM** currently interfaces to the following libraries:


.. hlist::
    :columns: 3
	 
    * `ARPACK <https://www.caam.rice.edu/software/ARPACK/>`_
    * `BLAS <http://www.netlib.org/blas/>`_
    * `OpenBLAS <http://www.openblas.net/>`_
    * `FFTW 3.3.2 <http://www.fftw.org>`_
    * `Ipopt 3.12.4 <https://github.com/coin-or/Ipopt>`_
    * `Gmm++ 4.2 <http://getfem.org/gmm.html>`_
    * `freeYams <https://www.ljll.math.upmc.fr/frey/software.html>`_
    * `METIS <http://glaros.dtc.umn.edu/gkhome/metis/metis/overview>`_
    * `ParMETIS <http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview>`_
    * `Mmg <https://www.mmgtools.org/>`_
    * `mshmet <https://www.ljll.math.upmc.fr/frey/software.html>`_
    * `MUMPS <http://mumps.enseeiht.fr/>`_
    * `NLopt 2.2.4 <http://ab-initio.mit.edu/wiki/index.php/NLopt>`_
    * `ScaLAPACK <http://www.netlib.org/scalapack/>`_
    * `Scotch <https://gforge.inria.fr/projects/scotch/>`_
    * `SuiteSparse <http://faculty.cse.tamu.edu/davis/suitesparse.html>`_
    * `SuperLU <http://crd-legacy.lbl.gov/~xiaoye/SuperLU/>`_
    * `TetGen <http://www.tetgen.org/>`_
    * `PETSc <https://www.mcs.anl.gov/petsc/>`_
    * `SLEPc <https://slepc.upv.es/>`_
    * `HTool <https://pierremarchand.netlify.com/project/htool/>`_
    * `HPDDM <https://github.com/hpddm/hpddm>`_
    * `BemTool <https://github.com/PierreMarchand20/BemTool>`_
    * `ParMmg <https://github.com/MmgTools/ParMmg>`_
	

|
|


.. image:: ../_static/img/Logo.png
   :align: center
   :width: 65%
   :alt: FreeFem++

.. _FreeFem++: https://freefem.org


|
|
|
|
|
|


.. only:: html
.. |logo1| image:: ../_static/img/logo_cnrs_SU_UP.png
   :scale: 40%
.. |logo2| image:: ../_static/img/inr_logo_rouge_300.jpg
   :scale: 20%
.. |logo3| image:: ../_static/img/logo_LJLL.png
   :scale: 60%
.. |logo4| image:: ../_static/img/logo_alpines.png
   :scale: 45%

.. centered::
  |logo1|  |logo2|  |logo3|    |logo4|
  



.. toctree::
   
   new-features
   installation
   download
   history
   citation
   authors
   contributing
