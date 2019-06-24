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
  |logo1|  |logo2|  |logo3|  |logo4|
  
|
|
|
|

.. raw:: html

    <style> .blue {color:blue;font-weight : bold;font-size:35px} </style>

.. role:: blue

:blue:`The FreeFEM` |version| :blue:`is out!`
	
.. only:: html

  .. raw:: html

    <div id="PDFbanner">
      <p>
        The user manuel is available <a href="/pdf/FreeFEM-documentation.pdf"><img src="../_static/img/logo_pdf.png"width="15" height=""></a>
      </p>
	  <p>
	    <b>Professor Frédéric Hecht</b>'s personal <a href="https://www.ljll.math.upmc.fr/hecht/" target="_blank">website</a> <br> 
	    
	    <font size="2">Laboratoire Jacques-Louis Lions, Université Pierre et Marie Curie, Paris</font>
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
	
|
|
|
|

Introduction
============

**FreeFEM** is a partial differential equation solver for non-linear multi-physics systems in 2D and 3D.

Problems involving partial differential equations from several branches of physics, such as fluid-structure interactions, require interpolations of data on several meshes and their manipulation within one program.
**FreeFEM** includes a fast interpolation algorithm and a language for the manipulation of data on multiple meshes.

**FreeFEM** is written in C++ and its language is a C++ idiom.

.. image:: ../_static/img/Logo.png
   :align: center
   :width: 65%
   :alt: FreeFem++

.. _FreeFem++: https://freefem.org

.. toctree::

   citation
   history
   new-features
   download
   installation
   authors
   contributing
