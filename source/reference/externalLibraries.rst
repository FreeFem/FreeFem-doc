.. role:: freefem(code)
  :language: freefem

.. _externalLibraries:

External libraries
==================

aniso
-----

boundaniso
~~~~~~~~~~

.. todo:: todo

BEC
---

BECtrap
~~~~~~~

.. todo:: todo

GPvortex
~~~~~~~~

.. todo:: todo

dxGPVortex
~~~~~~~~~~

.. todo:: todo

dyGPVortex
~~~~~~~~~~

.. todo:: todo

Binary I/O
----------

LoadVec
~~~~~~~

.. todo:: todo

LoadFlag
~~~~~~~~

.. todo:: todo

SaveVec
~~~~~~~

.. todo:: todo

flag
~~~~

.. todo:: todo

buildlayer
----------

buildlayers
~~~~~~~~~~~

.. todo:: todo

ClosePoints
-----------

radiusSearch
~~~~~~~~~~~~

.. todo:: todo

Voisinage
~~~~~~~~~

.. todo:: todo

neighborhood
~~~~~~~~~~~~

.. todo:: todo

ClosePoints2
~~~~~~~~~~~~

.. todo:: todo

ClosePoint
~~~~~~~~~~

.. todo:: todo

ClosePoints1
~~~~~~~~~~~~

.. todo:: todo

Curvature
---------

extractborder
~~~~~~~~~~~~~

Extract a border of a mesh.

.. code-block:: freefem
   :linenos:

   int Res = extractborder(Th, Label, Points);

Parameters:

-  ``Th`` (:freefem:`mesh` or :freefem:`mesh3`)
-  ``Label`` (:freefem:`int`) Label of the border to extract
-  ``Points`` (:freefem:`real[int, int]`) Extracted points Must be allocated as :freefem:`real[int, int] Points(3, 1);`

Output:

-  ``Res`` (:freefem:`real`) Length of the extracted border

curvature
~~~~~~~~~

.. todo:: todo

raxicurvature
~~~~~~~~~~~~~

.. todo:: todo

curves
~~~~~~

.. todo:: todo

setecurveabcisse
~~~~~~~~~~~~~~~~

.. todo:: todo

equiparameter
~~~~~~~~~~~~~

.. todo:: todo

Tresca
~~~~~~

.. todo:: todo

VonMises
~~~~~~~~

.. todo:: todo

dfft
----

Refer to the `FFTW documentation <http://www.fftw.org/>`__ for more informations.

plandfft
~~~~~~~~

.. todo:: todo

execute
~~~~~~~

.. todo:: todo

delete
~~~~~~

.. todo:: todo

dfft
~~~~

.. todo:: todo

map
~~~

.. todo:: todo

distance
--------

Need

.. code-block:: freefem
   :linenos:

   load "distance"

distance
~~~~~~~~

.. code-block:: freefem
   :linenos:

   distance(Th, d, dist, [distmax=DistMax]);

Parameters:

-  ``Th`` (:freefem:`mesh`)
-  ``d``
-  ``dist`` (:freefem:`real[int]`)

Output:

-

.. todo:: todo

checkdist
~~~~~~~~~

.. todo:: todo

DxWriter
--------

Dxaddmesh
~~~~~~~~~

.. todo:: todo

Dxaddtimeseries
~~~~~~~~~~~~~~~

.. todo:: todo

Dxaddsol2ts
~~~~~~~~~~~

.. todo:: todo

Element_P1bl
------------

expert
~~~~~~

.. todo:: todo

exactpartition
--------------

exactpartition
~~~~~~~~~~~~~~

.. todo:: todo

ff-AiryBiry
-----------

airy
~~~~

.. todo:: todo

biry
~~~~

.. todo:: todo

ff-cmaes
--------

cmaes
~~~~~

.. todo:: todo

.. _referenceFFGSLAWK:

ff_gsl_awk
----------

Refer to the `GSL documentation <https://www.gnu.org/software/gsl/doc/html/index.html>`__ for more informations

gslcdfugaussianP
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_ugaussian_P(a)

gslcdfugaussianQ
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_ugaussian_Q(a)

gslcdfugaussianPinv
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_ugaussian_Pinv(a)

gslcdfugaussianQinv
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_ugaussian_Qinv(a)

gslcdfgaussianP
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_gaussian_P(a, b)

gslcdfgaussianQ
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_gaussian_Q(a, b)

gslcdfgaussianPinv
~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_gaussian_Pinv(a, b)

gslcdfgaussianQinv
~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_gaussian_Qinv(a, b)

gslcdfgammaP
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_gamma_P(a, b, c)

gslcdfgammaQ
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_gamma_Q(a, b, c)

gslcdfgammaPinv
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_gamma_Pinv(a, b, c)

gslcdfgammaQinv
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_gamma_Pinv(a, b, c)

gslcdfcauchyP
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_cauchy_P(a, b)

gslcdfcauchyQ
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_cauchy_Q(a, b)

gslcdfcauchyPinv
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_cauchy_Pinv(a, b)

gslcdfcauchyQinv
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_cauchy_Qinv(a, b)

gslcdflaplaceP
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_lapalce_P(a, b)

gslcdflaplaceQ
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_lapalce_Q(a, b)

gslcdflaplacePinv
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_lapalce_Pinv(a, b)

gslcdflaplaceQinv
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_lapalce_Qinv(a, b)

gslcdfrayleighP
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_rayleigh_P(a, b)

gslcdfrayleighQ
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_rayleigh_Q(a, b)

gslcdfrayleighPinv
~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_rayleigh_Pinv(a, b)

gslcdfrayleighQinv
~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_rayleigh_Qinv(a, b)

gslcdfchisqP
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_chisq_P(a, b)

gslcdfchisqQ
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_chisq_Q(a, b)

gslcdfchisqPinv
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_chisq_Pinv(a, b)

gslcdfchisqQinv
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_chisq_Qinv(a, b)

gslcdfexponentialP
~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_exponential_P(a, b)

gslcdfexponentialQ
~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_exponential_Q(a, b)

gslcdfexponentialPinv
~~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_exponential_Pinv(a, b)

gslcdfexponentialQinv
~~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_exponential_Qinv(a, b)

gslcdfexppowP
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_exppow_P(a, b, c)

gslcdfexppowQ
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_exppow_Q(a, b, c)

gslcdftdistP
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_t_dist_P(a, b)

gslcdftdistQ
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_t_dist_Q(a, b)

gslcdftdistPinv
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_t_dist_Pinv(a, b)

gslcdftdistQinv
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_t_dist_Qinv(a, b)

gslcdffdistP
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_fdist_P(a, b, c)

gslcdffdistQ
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_fdist_Q(a, b, c)

gslcdffdistPinv
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_fdist_Pinv(a, b, c)

gslcdffdistQinv
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_fdist_Qinv(a, b, c)

gslcdfbetaP
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_beta_P(a, b, c)

gslcdfbetaQ
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_beta_Q(a, b, c)

gslcdfbetaPinv
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_beta_Pinv(a, b, c)

gslcdfbetaQinv
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_beta_Qinv(a, b, c)

gslcdfflatP
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_flat_P(a, b, c)

gslcdfflatQ
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_flat_Q(a, b, c)

gslcdfflatPinv
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_flat_Pinv(a, b, c)

gslcdfflatQinv
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_flat_Qinv(a, b, c)

gslcdflognormalP
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_lognormal_P(a, b, c)

gslcdflognormalQ
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_lognormal_Q(a, b, c)

gslcdflognormalPinv
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_lognormal_Pinv(a, b, c)

gslcdflognormalQinv
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_lognormal_Qinv(a, b, c)

gslcdfgumbel1P
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_gumbel1_P(a, b, c)

gslcdfgumbel1Q
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_gumbel1_Q(a, b, c)

gslcdfgumbel1Pinv
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_gumbel1_Pinv(a, b, c)

gslcdfgumbel1Qinv
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_gumbel1_Qinv(a, b, c)

gslcdfgumbel2P
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_gumbel2_P(a, b, c)

gslcdfgumbel2Q
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_gumbel2_Q(a, b, c)

gslcdfgumbel2Pinv
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_gumbel2_Pinv(a, b, c)

gslcdfgumbel2Qinv
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_gumbel2_Qinv(a, b, c)

gslcdfweibullP
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_weibull_P(a, b, c)

gslcdfweibullQ
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_weibull_Q(a, b, c)

gslcdfweibullPinv
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_weibull_Pinv(a, b, c)

gslcdfweibullQinv
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_weibull_Qinv(a, b, c)

gslcdfparetoP
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_pareto_P(a, b, c)

gslcdfparetoQ
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_pareto_Q(a, b, c)

gslcdfparetoPinv
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_pareto_Pinv(a, b, c)

gslcdfparetoQinv
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_pareto_Qinv(a, b, c)

gslcdflogisticP
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_logistic_P(a, b)

gslcdflogisticQ
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_logistic_Q(a, b)

gslcdflogisticPinv
~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_logistic_Pinv(a, b)

gslcdflogisticQinv
~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_logistic_Qinv(a, b)

gslcdfbinomialP
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_binomial_P(a, b, c)

gslcdfbinomialQ
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_binomial_Q(a, b, c)

gslcdfpoissonP
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_poisson_P(a, b)

gslcdfpoissonQ
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_poisson_Q(a, b)

gslcdfgeometricP
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_geometric_P(a, b)

gslcdfgeometricQ
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_geometric_Q(a, b)

gslcdfnegativebinomialP
~~~~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_negative_binomial_P(a, b, c)

gslcdfnegativebinomialQ
~~~~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_negative_binomial_Q(a, b, c)

gslcdfpascalP
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_pascal_P(a, b, c)

gslcdfpascalQ
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_cdf_pascal_Q(a, b, c)

gslranbernoullipdf
~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_bernoulli_pdf(a, b)

gslranbeta
~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_beta(a, b, c)

gslranbetapdf
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_beta_pdf(a, b, c)

gslranbinomialpdf
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_binomial_pdf(a, b, c)

gslranexponential
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_exponential(a, b)

gslranexponentialpdf
~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_exponential_pdf(a, b)

gslranexppow
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_exppow(a, b, c)

gslranexppowpdf
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_exppow_pdf(a, b, c)

gslrancauchy
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_cauchy(a, b)

gslrancauchypdf
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_cauchy_pdf(a, b)

gslranchisq
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_chisq(a, b)

gslranchisqpdf
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_chisq_pdf(a, b)

gslranerlang
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_erlang(a, b, c)

gslranerlangpdf
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_erlang_pdf(a, b, c)

gslranfdist
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_fdist(a, b, c)

gslranfdistpdf
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_fdist_pdf(a, b, c)

gslranflat
~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_flat(a, b, c)

gslranflatpdf
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_flat_pdf(a, b, c)

gslrangamma
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_gamma(a, b, c)

gslrangammaint
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_gamma_int(a, b, c)

gslrangammapdf
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_gamma_pdf(a, b, c)

gslrangammamt
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_gamma_mt(a, b, c)

gslrangammaknuth
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_gamma_knuth(a, b, c)

gslrangaussian
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_gaussian(a, b)

gslrangaussianratiomethod
~~~~~~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_gaussian_ratio_method(a, b)

gslrangaussianziggurat
~~~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_gaussian_ziggurat(a, b)

gslrangaussianpdf
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_gaussian_pdf(a, b)

gslranugaussian
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_ugaussian(a)

gslranugaussianratiomethod
~~~~~~~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_ugaussian_ratio_method(a)

gslranugaussianpdf
~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_ugaussian_pdf(a)

gslrangaussiantail
~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_gaussian_tail(a, b, c)

gslrangaussiantailpdf
~~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_gaussian_tail_pdf(a, b, c)

gslranugaussiantail
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_ugaussian_tail(a, b)

gslranugaussiantailpdf
~~~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_ugaussian_tail_pdf(a, b)

gslranlandau
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_landau(a)

gslranlandaupdf
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_landau_pdf(a)

gslrangeometricpdf
~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_geometric_pdf(a, b)

gslrangumbel1
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_gumbel1(a, b, c)

gslrangumbel1pdf
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_gumbel1_pdf(a, b, c)

gslrangumbel2
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_gumbel2(a, b, c)

gslrangumbel2pdf
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_gumbel2_pdf(a, b, c)

gslranlogistic
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_logistic(a, b)

gslranlogisticpdf
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_logistic_pdf(a, b)

gslranlognormal
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_lognormal(a, b, c)

gslranlognormalpdf
~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_lognormal_pdf(a, b, c)

gslranlogarithmicpdf
~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_logarithmic_pdf(a, b)

gslrannegativebinomialpdf
~~~~~~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_negative_binomial_pdf(a, b, c)

gslranpascalpdf
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_pascal_pdf(a, b, c)

gslranpareto
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_pareto(a, b, c)

gslranparetopdf
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_pareto_pdf(a, b, c)

gslranpoissonpdf
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_poisson_pdf(a, b)

gslranrayleigh
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_rayleigh(a, b)

gslranrayleighpdf
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_rayleigh_pdf(a, b)

gslranrayleightail
~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_rayleigh_tail(a, b, c)

gslranrayleightailpdf
~~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_rayleigh_tail_pdf(a, b, c)

gslrantdist
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_tdsit(a, b)

gslrantdistpdf
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_tdsit_pdf(a, b)

gslranlaplace
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_laplace(a, b)

gslranlaplacepdf
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_laplace_pdf(a, b)

gslranlevy
~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_levy(a, b, c)

gslranweibull
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_weibull(a, b, c)

gslranweibullpdf
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_ran_weibull_pdf(a, b, c)

gslsfairyAi
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_airy_Ai(a, b)

gslsfairyBi
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_airy_Bi(a, b)

gslsfairyAiscaled
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_airy_Ai_scaled(a, b)

gslsfairyBiscaled
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_airy_Bi_scaled(a, b)

gslsfairyAideriv
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_airy_Ai_deriv(a, b)

gslsfairyBideriv
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_airy_Bi_deriv(a, b)

gslsfairyAiderivscaled
~~~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_airy_Ai_deriv_scaled(a, b)

gslsfairyBiderivscaled
~~~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_airy_Bi_deriv_scaled(a, b)

gslsfairyzeroAi
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_airy_Ai(a, b)

gslsfairyzeroBi
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_airy_aero_Bi(a)

gslsfairyzeroAideriv
~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_airy_aero_Ai_deriv(a)

gslsfairyzeroBideriv
~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_airy_aero_Bi_deriv(a)

gslsfbesselJ0
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_J0(a)

gslsfbesselJ1
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_J1(a)

gslsfbesselJn
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_Jn(a, b)

gslsfbesselY0
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_Y0(a)

gslsfbesselY1
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_Y1(a)

gslsfbesselYn
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_Yn(a, b)

gslsfbesselI0
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_I0(a)

gslsfbesselI1
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_I1(a)

gslsfbesselIn
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_In(a, b)

gslsfbesselI0scaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_I0_scaled(a)

gslsfbesselI1scaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_I1_scaled(a)

gslsfbesselInscaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_In_scaled(a, b)

gslsfbesselK0
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_K0(a)

gslsfbesselK1
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_K1(a)

gslsfbesselKn
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_Kn(a, b)

gslsfbesselK0scaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_K0_scaled(a)

gslsfbesselK1scaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_K1_scaled(a)

gslsfbesselKnscaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_Kn_scaled(a, b)

gslsfbesselj0
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_j0(a)

gslsfbesselj1
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_j1(a)

gslsfbesselj2
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_j2(a)

gslsfbesseljl
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_jl(a, b)

gslsfbessely0
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_y0(a)

gslsfbessely1
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_y0(a)

gslsfbessely2
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_y0(a)

gslsfbesselyl
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_jl(a, b)

gslsfbesseli0scaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_i0_scaled(a)

gslsfbesseli1scaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_i1_scaled(a)

gslsfbesseli2scaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_i2_scaled(a)

gslsfbesselilscaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_il_scaled(a, b)

gslsfbesselk0scaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_k0_scaled(a)

gslsfbesselk1scaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_k1_scaled(a)

gslsfbesselk2scaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_k2_scaled(a)

gslsfbesselklscaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_kl_scaled(a, b)

gslsfbesselJnu
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_Jnu(a, b)

gslsfbesselYnu
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_Ynu(a, b)

gslsfbesselInuscaled
~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_Inu_scaled(a, b)

gslsfbesselInu
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_Inu(a, b)

gslsfbesselKnuscaled
~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_Knu_scaled(a, b)

gslsfbesselKnu
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_Knu(a, b)

gslsfbessellnKnu
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_lnKnu(a, b)

gslsfbesselzeroJ0
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_zero_J0(a)

gslsfbesselzeroJ1
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_zero_J1(a)

gslsfbesselzeroJnu
~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_bessel_zero_Jnu(a, b)

gslsfclausen
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_clausen(a)

gslsfhydrogenicR1
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_hydrogenicR_1(a, b)

gslsfdawson
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_dawnson(a)

gslsfdebye1
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_debye_1(a)

gslsfdebye2
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_debye_2(a)

gslsfdebye3
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_debye_3(a)

gslsfdebye4
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_debye_4(a)

gslsfdebye5
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_debye_5(a)

gslsfdebye6
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_debye_6(a)

gslsfdilog
~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_dilog(a)

gslsfmultiply
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_multiply(a, b)

gslsfellintKcomp
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_ellint_Kcomp(a, b)

gslsfellintEcomp
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_ellint_Ecomp(a, b)

gslsfellintPcomp
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_ellint_Pcomp(a, b, c)

gslsfellintDcomp
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_ellint_Dcomp(a, b)

gslsfellintF
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_ellint_F(a, b, c)

gslsfellintE
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_ellint_E(a, b, c)

gslsfellintRC
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_ellint_RC(a, b, c)

gslsferfc
~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_erfc(a)

gslsflogerfc
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_log_erfc(a)

gslsferf
~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_erf(a)

gslsferfZ
~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_erf_Z(a)

gslsferfQ
~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_erf_Q(a)

gslsfhazard
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_hazard(a)

gslsfexp
~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_exp(a)

gslsfexpmult
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_exp_mult(a, b)

gslsfexpm1
~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_expm1(a)

gslsfexprel
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_exprel(a)

gslsfexprel2
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_exprel_2(a)

gslsfexpreln
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_exprel_n(a, b)

gslsfexpintE1
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_expint_E1(a)

gslsfexpintE2
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_expint_E2(a)

gslsfexpintEn
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_expint_En(a, b)

gslsfexpintE1scaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_expint_E1_scaled(a)

gslsfexpintE2scaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_expint_E1_scaled(a)

gslsfexpintEnscaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_expint_En_scaled(a, b)

gslsfexpintEi
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_expint_Ei(a)

gslsfexpintEiscaled
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_expint_Ei_scaled(a)

gslsfShi
~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_Shi(a)

gslsfChi
~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_Chi(a)

gslsfexpint3
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_expint_3(a)

gslsfSi
~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_Si(a)

gslsfCi
~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_Ci(a)

gslsfatanint
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_atanint(a)

gslsffermidiracm1
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_fermi_dirac_m1(a)

gslsffermidirac0
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_fermi_dirac_0(a)

gslsffermidirac1
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_fermi_dirac_1(a)

gslsffermidirac2
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_fermi_dirac_2(a)

gslsffermidiracint
~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_fermi_dirac_int(a, b)

gslsffermidiracmhalf
~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_fermi_dirac_mhalf(a)

gslsffermidirachalf
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_fermi_dirac_half(a)

gslsffermidirac3half
~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_fermi_dirac_3half(a)

gslsffermidiracinc0
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_fermi_dirac_inc_0(a, b)

gslsflngamma
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_lngamma(a)

gslsfgamma
~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_gamma(a)

gslsfgammastar
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_gammastar(a)

gslsfgammainv
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_gammainv(a)

gslsftaylorcoeff
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_taylorcoeff(a, b)

gslsffact
~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_fact(a)

gslsfdoublefact
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_doublefact(a)

gslsflnfact
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_lnfact(a)

gslsflndoublefact
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_lndoublefact(a)

gslsflnchoose
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_lnchoose(a, b)

gslsfchoose
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_choose(a, b)

gslsflnpoch
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_lnpoch(a, b)

gslsfpoch
~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_poch(a, b)

gslsfpochrel
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_pochrel(a, b)

gslsfgammaincQ
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_gamma_inc_Q(a, b)

gslsfgammaincP
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_gamma_inc_P(a, b)

gslsfgammainc
~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_gamma_inc(a, b)

gslsflnbeta
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_lnbeta(a, b)

gslsfbeta
~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_beta(a, b)

gslsfbetainc
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_betaçinc(a, b, c)

gslsfgegenpoly1
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_gegenpoly_1(a, b)

gslsfgegenpoly2
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_gegenpoly_2(a, b)

gslsfgegenpoly3
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_gegenpoly_3(a, b)

gslsfgegenpolyn
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_gegenpoly_n(a, b, c)

gslsfhyperg0F1
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_hyperg_0F1(a, b)

gslsfhyperg1F1int
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_hyperg_1F1_inc(a, b, c)

gslsfhyperg1F1
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_hyperg_1F1(a, b, c)

gslsfhypergUint
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_hyperg_U_inc(a, b, c)

gslsfhypergU
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_hyperg_U(a, b, c)

gslsfhyperg2F0
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_hyperg_U_2F0(a, b, c)

gslsflaguerre1
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_laguerre_1(a, b)

gslsflaguerre2
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_laguerre_2(a, b)

gslsflaguerre3
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_laguerre_3(a, b)

gslsflaguerren
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_laguerre_n(a, b, c)

gslsflambertW0
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_lambert_W0(a)

gslsflambertWm1
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_lambert_Wm1(a)

gslsflegendrePl
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_legendre_Pl(a, b)

gslsflegendreP1
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_legendre_P1(a)

gslsflegendreP2
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_legendre_P2(a)

gslsflegendreP3
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_legendre_P3(a)

gslsflegendreQ0
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_legendre_Q0(a)

gslsflegendreQ1
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_legendre_Q1(a)

gslsflegendreQl
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_legendre_Ql(a, b)

gslsflegendrePlm
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_legendre_Plm(a, b, c)

gslsflegendresphPlm
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_legendre_sphP1m(a, b, c)

gslsflegendrearraysize
~~~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_legendre_array_size(a, b)

gslsfconicalPhalf
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_conicalP_half(a, b)

gslsfconicalPmhalf
~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_conicalP_mhalf(a, b)

gslsfconicalP0
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_conicalP_0(a, b)

gslsfconicalP1
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_conicalP_1(a, b)

gslsfconicalPsphreg
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_conicalP_sph_reg(a, b, c)

gslsfconicalPcylreg
~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_conicalP_cyl_reg(a, b, c)

gslsflegendreH3d0
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_legendre_H3d_0(a, b)

gslsflegendreH3d1
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_legendre_H3d_1(a, b)

gslsflegendreH3d
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_legendre_H3d(a, b, c)

gslsflog
~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_log(a)

gslsflogabs
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_log_abs(a)

gslsflog1plusx
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_log_1plusx(a)

gslsflog1plusxmx
~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_log_1plusx_mx(a)

gslsfpowint
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_pow_int(a, b)

gslsfpsiint
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_psi_int(a)

gslsfpsi
~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_psi(a)

gslsfpsi1piy
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_psi_1piy(a)

gslsfpsi1int
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_psi_1_int(a)

gslsfpsi1
~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_psi_1(a)

gslsfpsin
~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_psi_n(a, b)

gslsfsynchrotron1
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_synchrotron_1(a)

gslsfsynchrotron2
~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_synchrotron_2(a)

gslsftransport2
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_transport_2(a)

gslsftransport3
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_transport_3(a)

gslsftransport4
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_transport_4(a)

gslsftransport5
~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_transport_5(a)

gslsfsin
~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_sin(a)

gslsfcos
~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_cos(a)

gslsfhypot
~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_hypot(a, b)

gslsfsinc
~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_sinc(a)

gslsflnsinh
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_lnsinh(a)

gslsflncosh
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_lncosh(a)

gslsfanglerestrictsymm
~~~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_andle_restrict_symm(a)

gslsfanglerestrictpos
~~~~~~~~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_angle_restrict_pos(a)

gslsfzetaint
~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_zeta_int(a)

gslsfzeta
~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_zeta(a)

gslsfzetam1
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_zetam1(a)

gslsfzetam1int
~~~~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_zetam1_int(a)

gslsfhzeta
~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_hzeta(a, b)

gslsfetaint
~~~~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_eta_int(a)

gslsfeta
~~~~~~~~

Link to:

.. code-block:: cpp
   :linenos:

   gsl_sf_eta(a)

ff-Ipopt
--------

Refer to the `Ipopt documentation <https://projects.coin-or.org/Ipopt>`__ for more informations.

IPOPT
~~~~~

.. todo:: todo

fflapack
--------

Refer to the `LAPACK documentation <http://www.netlib.org/lapack/>`__ for more informations.

inv
~~~

.. todo:: todo

dgeev
~~~~~

.. todo:: todo

zgeev
~~~~~

.. todo:: todo

geev
~~~~

.. todo:: todo

geev
~~~~

.. todo:: todo

dggev
~~~~~

.. todo:: todo

zggev
~~~~~

.. todo:: todo

dsygvd
~~~~~~

.. todo:: todo

dgesdd
~~~~~~

.. todo:: todo

zhegv
~~~~~

.. todo:: todo

dsyev
~~~~~

.. todo:: todo

zheev
~~~~~

.. todo:: todo

ff-mmap-semaphore
-----------------

Wait
~~~~

.. todo:: todo

trywait
~~~~~~~

.. todo:: todo

Post
~~~~

.. todo:: todo

msync
~~~~~

.. todo:: todo

Read
~~~~

.. todo:: todo

Write
~~~~~

.. todo:: todo

ffnewuoa
--------

newuoa
~~~~~~

.. todo:: todo

ff-NLopt
--------

Refer to the `NLOPT documentation <https://nlopt.readthedocs.io/en/latest/>`__ for more informations.

nloptDIRECT
~~~~~~~~~~~

.. todo:: todo

nloptDIRECTL
~~~~~~~~~~~~

.. todo:: todo

nloptDIRECTLRand
~~~~~~~~~~~~~~~~

.. todo:: todo

nloptDIRECTScal
~~~~~~~~~~~~~~~

.. todo:: todo

nloptDIRECTNoScal
~~~~~~~~~~~~~~~~~

.. todo:: todo

nloptDIRECTLNoScal
~~~~~~~~~~~~~~~~~~

.. todo:: todo

nloptDIRECTLRandNoScal
~~~~~~~~~~~~~~~~~~~~~~

.. todo:: todo

nloptOrigDIRECT
~~~~~~~~~~~~~~~

.. todo:: todo

nloptOrigDIRECTL
~~~~~~~~~~~~~~~~

.. todo:: todo

nloptStoGO
~~~~~~~~~~

.. todo:: todo

nloptStoGORand
~~~~~~~~~~~~~~

.. todo:: todo

nloptLBFGS
~~~~~~~~~~

.. todo:: todo

nloptPRAXIS
~~~~~~~~~~~

.. todo:: todo

nloptVar1
~~~~~~~~~

.. todo:: todo

nloptVar2
~~~~~~~~~

.. todo:: todo

nloptTNewton
~~~~~~~~~~~~

.. todo:: todo

nloptTNewtonRestart
~~~~~~~~~~~~~~~~~~~

.. todo:: todo

nloptTNewtonPrecond
~~~~~~~~~~~~~~~~~~~

.. todo:: todo

nloptTNewtonPrecondRestart
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. todo:: todo

nloptCRS2
~~~~~~~~~

.. todo:: todo

nloptMMA
~~~~~~~~

.. todo:: todo

nloptCOBYLA
~~~~~~~~~~~

.. todo:: todo

nloptNEWUOA
~~~~~~~~~~~

.. todo:: todo

nloptNEWUOABound
~~~~~~~~~~~~~~~~

.. todo:: todo

nloptNelderMead
~~~~~~~~~~~~~~~

.. todo:: todo

nloptSbplx
~~~~~~~~~~

.. todo:: todo

nloptBOBYQA
~~~~~~~~~~~

.. todo:: todo

nloptISRES
~~~~~~~~~~

.. todo:: todo

nloptSLSQP
~~~~~~~~~~

.. todo:: todo

nloptMLSL
~~~~~~~~~

.. todo:: todo

nloptMLSLLDS
~~~~~~~~~~~~

.. todo:: todo

nloptAUGLAG
~~~~~~~~~~~

.. todo:: todo

nloptAUGLAGEQ
~~~~~~~~~~~~~

.. todo:: todo

ffrandom
--------

srandomdev
~~~~~~~~~~

.. todo:: todo

srandom
~~~~~~~

.. todo:: todo

random
~~~~~~

.. todo:: todo

FreeFemQA
---------

MeshGenQA
~~~~~~~~~

.. todo:: todo

freeyams
--------

freeyams
~~~~~~~~

.. todo:: todo


gmsh
----

Need

.. code-block:: freefem
   :linenos:

   load "gsmh"

The ``gmsh`` software is available `here <http://gmsh.info/>`__

gmshload
~~~~~~~~

Load a 2D mesh build with Gmsh.

.. code-block:: freefem
   :linenos:

   mesh Th = gmshload(MeshFile, [reftri=RefTri], [renum=Renum]);

Parameters:

-  ``MeshFile`` (:freefem:`string`) Mesh file name
-  :freefem:`reftri=` (.. todo:: todo)
-  :freefem:`renum=` (.. todo:: todo)

Output:

-  ``Th`` (:freefem:`mesh`)

gmshload3
~~~~~~~~~

Load a 3D mesh build with Gmsh.

.. code-block:: freefem
   :linenos:

   mesh3 Th = gmshload(MeshFile, [reftet=RefTet], [renum=Renum]);

Parameters:

-  ``MeshFile`` (:freefem:`string`) Mesh file name
-  :freefem:`reftet=` (.. todo:: todo)
-  :freefem:`renum=` (.. todo:: todo)

Output:

-  ``Th`` (:freefem:`mesh3`)

savegmsh
~~~~~~~~

.. todo:: todo

gsl
---

gslpolysolvequadratic
~~~~~~~~~~~~~~~~~~~~~

.. todo:: todo

gslpolysolvecubic
~~~~~~~~~~~~~~~~~

.. todo:: todo

gslpolycomplexsolve
~~~~~~~~~~~~~~~~~~~

.. todo:: todo

gslrnguniform
~~~~~~~~~~~~~

.. todo:: todo

gslrnguniformpos
~~~~~~~~~~~~~~~~

.. todo:: todo

gslname
~~~~~~~

.. todo:: todo

gslrngget
~~~~~~~~~

.. todo:: todo

gslrngmin
~~~~~~~~~

.. todo:: todo

gslrngmax
~~~~~~~~~

.. todo:: todo

gslrngset
~~~~~~~~~

.. todo:: todo

gslrngtype
~~~~~~~~~~

.. todo:: todo

ilut
----

applyIlutPrecond
~~~~~~~~~~~~~~~~

.. todo:: todo

makeIlutPrecond
~~~~~~~~~~~~~~~

.. todo:: todo

iohdf5
------

savehdf5sol
~~~~~~~~~~~

.. todo:: todo

iovtk
-----

savevtk
~~~~~~~

Save mesh or solution in vtk/vtu format.

.. code-block:: freefem
   :linenos:

   savevtk(FileName, Th, [Ux, Uy, Uz], p, [dataname=DataName], [withsurfacemesh=WithSurfaceMesh], [order=Order], [floatmesh=FloatMesh], [floatsol=FloatSol], [bin=Bin], [swap=Swap]);

Parameters:

-  ``FileName`` (:freefem:`string`) File name: ``*.vtk`` or
   ``*.vtu``
-  ``Th`` (:freefem:`mesh` or :freefem:`mesh3`)
-  ``[Ux, Uy, Uz], p`` (:freefem:`fespace` function of vector of :freefem:`fespace` functions) Solutions to save, as much as wanted
-  :freefem:`dataname=` (:freefem:`string`) Name of solutions, seprated by a space
-  :freefem:`withsurfacemesh=` (:freefem:`bool`)
   .. todo:: todo
-  :freefem:`order=` (:freefem:`int[int]`) Order of solutions.

   Available: 0 or 1
-  :freefem:`floatmesh=` (:freefem:`bool`) .. todo:: todo
-  :freefem:`floatsol=` (:freefem:`bool`) .. todo:: todo
-  :freefem:`bin=` (:freefem:`bool`) If true, save file in binary format
-  :freefem:`swap` (:freefem:`bool`) .. todo:: todo

Output:

-  None

vtkload
~~~~~~~

.. todo:: todo

vtkload3
~~~~~~~~

.. todo:: todo

isoline
-------

Need

.. code-block:: freefem
   :linenos:

   load "isoline"

isoline
~~~~~~~

.. code-block:: freefem
   :linenos:

   int N = isoline(Th, u, xy, iso=Iso, close=Close, smoothing=Smoothing, ratio=Ratio, eps=Eps, beginend=BeginEnd, file=File);

.. todo:: todo

Curve
~~~~~

.. todo:: todo

Area
~~~~

.. todo:: todo

findallocalmin
~~~~~~~~~~~~~~

.. todo:: todo

lapack
------

inv
~~~

.. todo:: todo

dgeev
~~~~~

.. todo:: todo

zgeev
~~~~~

.. todo:: todo

geev
~~~~

.. todo:: todo

dggev
~~~~~

.. todo:: todo

zggev
~~~~~

.. todo:: todo

dsygvd
~~~~~~

.. todo:: todo

dgesdd
~~~~~~

.. todo:: todo

zhegv
~~~~~

.. todo:: todo

dsyev
~~~~~

.. todo:: todo

zheev
~~~~~

.. todo:: todo

dgelsy
~~~~~~

.. todo:: todo

lgbmo
-----

bmo
~~~

.. todo:: todo

mat_dervieux
------------

MatUpWind1
~~~~~~~~~~

.. todo:: todo

mat_psi
-------

MatUpWind0
~~~~~~~~~~

.. todo:: todo

.. _referenceMedit:

medit
-----

medit
~~~~~

.. todo:: todo

savesol
~~~~~~~

.. todo:: todo

readsol
~~~~~~~

.. todo:: todo

metis
-----

metisnodal
~~~~~~~~~~

.. todo:: todo

metisdual
~~~~~~~~~

.. todo:: todo

MetricKuate
-----------

MetricKuate
~~~~~~~~~~~

.. todo:: todo

MetricPk
--------

MetricPk
~~~~~~~~

.. todo:: todo

mmg3d
-----

mmg3d
~~~~~

.. todo:: todo

mmg3d-v4.0
----------

mmg3d
~~~~~

.. todo:: todo

msh3
----

change
~~~~~~

.. todo:: todo

movemesh23
~~~~~~~~~~

.. todo:: todo

movemesh2D3Dsurf
~~~~~~~~~~~~~~~~

.. todo:: todo

movemesh3
~~~~~~~~~

.. todo:: todo

movemesh
~~~~~~~~

.. todo:: todo

movemesh3D
~~~~~~~~~~

.. todo:: todo

deplacement
~~~~~~~~~~~

.. todo:: todo

checkbemesh
~~~~~~~~~~~

.. todo:: todo

buildlayers
~~~~~~~~~~~

.. todo:: todo

bcube
~~~~~

.. todo:: todo

cube
~~~~

Construct a cubic mesh.

.. code-block:: freefem
   :linenos:

   mesh3 Th = cube(nnX, nnY, nnZ, [X(x), Y(y), Z(z)], [label=Label], [flags=Flags], [region=Region]);

Parameters:

-  ``nnX`` (:freefem:`int`) Number of discretization point along :math:`x`
-  ``nnY`` (:freefem:`int`) Number of discretization point along :math:`y`
-  ``nnZ`` (:freefem:`int`) Number of discretization point along :math:`z`
-  ``X(x)`` (:freefem:`func`) *[Optional]*\  Affine function of :math:`x` to define the length Default: ``x``
-  ``Y(y)`` (:freefem:`func`) *[Optional]*\  Affine function of :math:`y` to define the width Default: ``y``
-  ``Z(z)`` (:freefem:`func`) *[Optional]*\  Affine function of :math:`z` to define the height Default: ``z``
-  :freefem:`label=` (:freefem:`int[int]`) *[Optional]*

   List of surface labels Default: ``[1, 2, 3, 4, 5, 6]``
-  :freefem:`flags=` (:freefem:`int`) *[Optional]*

   Refer to :ref:`square <functionSquare>`
-  :freefem:`region=` (:freefem:`int`) *[Optional]*

   Region number of the cube volume Default: ``0``

Output:

-  ``Th`` (:freefem:`mesh3`) Cube mesh

trunc
~~~~~

.. todo:: todo

gluemesh
~~~~~~~~

.. todo:: todo

extract
~~~~~~~

.. todo:: todo

showborder
~~~~~~~~~~

.. todo:: todo

getborder
~~~~~~~~~

.. todo:: todo

AddLayers
~~~~~~~~~

.. todo:: todo

mshmet
------

mshmet
~~~~~~

.. todo:: todo

MUMPS
-----

defaulttoMUMPSseq
~~~~~~~~~~~~~~~~~

.. todo:: todo

MUMPS_seq
---------

defaulttoMUMPSseq
~~~~~~~~~~~~~~~~~

.. todo:: todo



netgen
------

netg
~~~~

.. todo:: todo

netgstl
~~~~~~~

.. todo:: todo

netgload
~~~~~~~~

.. todo:: todo

NewSolver
---------

defaulttoUMFPACK
~~~~~~~~~~~~~~~~

.. todo:: todo

PARDISO
-------

defaulttoPARDISO
~~~~~~~~~~~~~~~~

.. todo:: todo

ompsetnumthreads
~~~~~~~~~~~~~~~~

.. todo:: todo

ompgetnumthreads
~~~~~~~~~~~~~~~~

.. todo:: todo

ompgetmaxthreads
~~~~~~~~~~~~~~~~

.. todo:: todo

pcm2rnm
-------

readpcm
~~~~~~~

.. todo:: todo

pipe
----

flush
~~~~~

.. todo:: todo

sleep
~~~~~

.. todo:: todo

usleep
~~~~~~

.. todo:: todo

qf11to25
--------

QF1d
~~~~

.. todo:: todo

QF2d
~~~~

.. todo:: todo

QF3d
~~~~

.. todo:: todo

tripleQF
~~~~~~~~

scotch
------

scotch
~~~~~~

.. todo:: todo

shell
-----

readdir
~~~~~~~

.. todo:: todo

unlink
~~~~~~

.. todo:: todo

rmdir
~~~~~

.. todo:: todo

cddir
~~~~~

.. todo:: todo

chdir
~~~~~

.. todo:: todo

basename
~~~~~~~~

.. todo:: todo

dirname
~~~~~~~

.. todo:: todo

mkdir
~~~~~

.. todo:: todo

chmod
~~~~~

.. todo:: todo

cpfile
~~~~~~

.. todo:: todo

stat
~~~~

.. todo:: todo

isdir
~~~~~

.. todo:: todo

getenv
~~~~~~

.. todo:: todo

setenv
~~~~~~

.. todo:: todo

unsetenv
~~~~~~~~

.. todo:: todo

splitedges
----------

SplitedgeMesh
~~~~~~~~~~~~~

.. todo:: todo

splitmesh12
-----------

splitmesh12
~~~~~~~~~~~

.. todo:: todo

splitmesh3
----------

splitmesh3
~~~~~~~~~~

.. todo:: todo

splitmesh4
----------

splimesh4
~~~~~~~~~

.. todo:: todo

splitmesh6
----------

splitmesh6
~~~~~~~~~~

.. todo:: todo

SuperLu
-------

defaulttoSuperLu
~~~~~~~~~~~~~~~~

.. todo:: todo

symmetrizeCSR
-------------

symmetrizeCSR
~~~~~~~~~~~~~

.. todo:: todo

tetgen
------

Refer to the `Tetgen documentation <http://wias-berlin.de/software/tetgen/>`__ for more informations.

tetgconvexhull
~~~~~~~~~~~~~~

.. todo:: todo

tetgtransfo
~~~~~~~~~~~

.. todo:: todo

tetg
~~~~

Build a 3D mesh from a surface.

.. code-block:: freefem
   :linenos:

   mesh3 Th = tetg(Th0, [reftet=RefTet], [label=Label], [switch=Switch], [nbofholes=NbOfHoles], [holelist=HoleList], [nbofregions=NbOfRegions], [regionlist=RegionList], [nboffacetcl=NbOfFaceTcl], [facetcl=FaceTcl])

.. todo:: todo

tetgreconstruction
~~~~~~~~~~~~~~~~~~

.. todo:: todo

UMFPACK64
---------

defaulttoUMFPACK64
~~~~~~~~~~~~~~~~~~

.. todo:: todo

VTK_writer_3d
-------------

Vtkaddmesh
~~~~~~~~~~

.. todo:: todo

Vtkaddscalar
~~~~~~~~~~~~

.. todo:: todo

VTK_writer
----------

Vtkaddmesh
~~~~~~~~~~

.. todo:: todo

Vtkaddscalar
~~~~~~~~~~~~
