<!--
## function

```freefem
example code
```

<u>Parameters:</u>

-

<u>Output:</u>

-
-->

## aniso

### boundaniso

$\codered$

## BEC

### BECtrap

$\codered$

### GPvortex

$\codered$

### dxGPVortex

$\codered$

### dyGPVortex

$\codered$

## Binary I/O

### LoadVec

$\codered$

### LoadFlag

$\codered$

### SaveVec

$\codered$

### flag

$\codered$

## buildlayer

### buildlayers

$\codered$

## ClosePoints

### radiusSearch

$\codered$

### Voisinage

$\codered$

### neighborhood

$\codered$

### ClosePoints2

$\codered$

### ClosePoint

$\codered$

### ClosePoints1

$\codered$

## Curvature

### extractborder
Extract a border of a mesh.

```freefem
int Res = extractborder(Th, Label, Points);
```

<u>Parameters:</u>

 - `Th` (`:::freefem mesh` or `:::freefem mesh3`)
 - `Label` (`:::freefem int`)<br/>
 Label of the border to extract
 - `Points` (`:::freefem real[int, int]`)<br/>
 Extracted points<br/>
 Must be allocated as `:::freefem real[int, int] Points(3, 1);`

<u>Output:</u>

 - `Res` (`:::freefem real`)<br/>
 Length of the extracted border

### curvature

$\codered$

### raxicurvature

$\codered$

### curves

$\codered$

### setecurveabcisse

$\codered$

### equiparameter

$\codered$

### Tresca

$\codered$

### VonMises

$\codered$

## dfft

Refer to the [FFTW documentation](http://www.fftw.org/) for more informations.

### plandfft

$\codered$

### execute

$\codered$

### delete

$\codered$

### dfft

$\codered$

### map

$\codered$

## distance

Need
```freefem
load "distance"
```

### distance

```freefem
distance(Th, d, dist, [distmax=DistMax]);
```

<u>Parameters:</u>

- `Th` (`:::freefem mesh`)
- `d`
- `dist` (`:::freefem real[int]`)

<u>Output:</u>

-

$\codered$

### checkdist

$\codered$

## DxWriter

### Dxaddmesh

$\codered$

### Dxaddtimeseries

$\codered$

### Dxaddsol2ts

$\codered$

## Element_P1bl

### expert

$\codered$

## exactpartition

### exactpartition

$\codered$

## ff-AiryBiry

### airy

$\codered$

### biry

$\codered$

## ff-cmaes

### cmaes

$\codered$

## ff_gsl_awk

Refer to the [GSL documentation](https://www.gnu.org/software/gsl/doc/html/index.html) for more informations

### gslcdfugaussianP
Link to:
```cpp
gsl_cdf_ugaussian_P(a)
```

### gslcdfugaussianQ
Link to:
```cpp
gsl_cdf_ugaussian_Q(a)
```

### gslcdfugaussianPinv
Link to:
```cpp
gsl_cdf_ugaussian_Pinv(a)
```

### gslcdfugaussianQinv
Link to:
```cpp
gsl_cdf_ugaussian_Qinv(a)
```

### gslcdfgaussianP
Link to:
```cpp
gsl_cdf_gaussian_P(a, b)
```

### gslcdfgaussianQ
Link to:
```cpp
gsl_cdf_gaussian_Q(a, b)
```

### gslcdfgaussianPinv
Link to:
```cpp
gsl_cdf_gaussian_Pinv(a, b)
```

### gslcdfgaussianQinv
Link to:
```cpp
gsl_cdf_gaussian_Qinv(a, b)
```

### gslcdfgammaP
Link to:
```cpp
gsl_cdf_gamma_P(a, b, c)
```

### gslcdfgammaQ
Link to:
```cpp
gsl_cdf_gamma_Q(a, b, c)
```

### gslcdfgammaPinv
Link to:
```cpp
gsl_cdf_gamma_Pinv(a, b, c)
```

### gslcdfgammaQinv
Link to:
```cpp
gsl_cdf_gamma_Pinv(a, b, c)
```

### gslcdfcauchyP
Link to:
```cpp
gsl_cdf_cauchy_P(a, b)
```

### gslcdfcauchyQ
Link to:
```cpp
gsl_cdf_cauchy_Q(a, b)
```

### gslcdfcauchyPinv
Link to:
```cpp
gsl_cdf_cauchy_Pinv(a, b)
```

### gslcdfcauchyQinv
Link to:
```cpp
gsl_cdf_cauchy_Qinv(a, b)
```

### gslcdflaplaceP
Link to:
```cpp
gsl_cdf_lapalce_P(a, b)
```

### gslcdflaplaceQ
Link to:
```cpp
gsl_cdf_lapalce_Q(a, b)
```

### gslcdflaplacePinv
Link to:
```cpp
gsl_cdf_lapalce_Pinv(a, b)
```

### gslcdflaplaceQinv
Link to:
```cpp
gsl_cdf_lapalce_Qinv(a, b)
```

### gslcdfrayleighP
Link to:
```cpp
gsl_cdf_rayleigh_P(a, b)
```

### gslcdfrayleighQ
Link to:
```cpp
gsl_cdf_rayleigh_Q(a, b)
```

### gslcdfrayleighPinv
Link to:
```cpp
gsl_cdf_rayleigh_Pinv(a, b)
```

### gslcdfrayleighQinv
Link to:
```cpp
gsl_cdf_rayleigh_Qinv(a, b)
```

### gslcdfchisqP
Link to:
```cpp
gsl_cdf_chisq_P(a, b)
```

### gslcdfchisqQ
Link to:
```cpp
gsl_cdf_chisq_Q(a, b)
```

### gslcdfchisqPinv
Link to:
```cpp
gsl_cdf_chisq_Pinv(a, b)
```

### gslcdfchisqQinv
Link to:
```cpp
gsl_cdf_chisq_Qinv(a, b)
```

### gslcdfexponentialP
Link to:
```cpp
gsl_cdf_exponential_P(a, b)
```

### gslcdfexponentialQ
Link to:
```cpp
gsl_cdf_exponential_Q(a, b)
```

### gslcdfexponentialPinv
Link to:
```cpp
gsl_cdf_exponential_Pinv(a, b)
```

### gslcdfexponentialQinv
Link to:
```cpp
gsl_cdf_exponential_Qinv(a, b)
```

### gslcdfexppowP
Link to:
```cpp
gsl_cdf_exppow_P(a, b, c)
```

### gslcdfexppowQ
Link to:
```cpp
gsl_cdf_exppow_Q(a, b, c)
```

### gslcdftdistP
Link to:
```cpp
gsl_cdf_t_dist_P(a, b)
```

### gslcdftdistQ
Link to:
```cpp
gsl_cdf_t_dist_Q(a, b)
```

### gslcdftdistPinv
Link to:
```cpp
gsl_cdf_t_dist_Pinv(a, b)
```

### gslcdftdistQinv
Link to:
```cpp
gsl_cdf_t_dist_Qinv(a, b)
```

### gslcdffdistP
Link to:
```cpp
gsl_cdf_fdist_P(a, b, c)
```

### gslcdffdistQ
Link to:
```cpp
gsl_cdf_fdist_Q(a, b, c)
```

### gslcdffdistPinv
Link to:
```cpp
gsl_cdf_fdist_Pinv(a, b, c)
```

### gslcdffdistQinv
Link to:
```cpp
gsl_cdf_fdist_Qinv(a, b, c)
```

### gslcdfbetaP
Link to:
```cpp
gsl_cdf_beta_P(a, b, c)
```

### gslcdfbetaQ
Link to:
```cpp
gsl_cdf_beta_Q(a, b, c)
```

### gslcdfbetaPinv
Link to:
```cpp
gsl_cdf_beta_Pinv(a, b, c)
```

### gslcdfbetaQinv
Link to:
```cpp
gsl_cdf_beta_Qinv(a, b, c)
```

### gslcdfflatP
Link to:
```cpp
gsl_cdf_flat_P(a, b, c)
```

### gslcdfflatQ
Link to:
```cpp
gsl_cdf_flat_Q(a, b, c)
```

### gslcdfflatPinv
Link to:
```cpp
gsl_cdf_flat_Pinv(a, b, c)
```

### gslcdfflatQinv
Link to:
```cpp
gsl_cdf_flat_Qinv(a, b, c)
```

### gslcdflognormalP
Link to:
```cpp
gsl_cdf_lognormal_P(a, b, c)
```

### gslcdflognormalQ
Link to:
```cpp
gsl_cdf_lognormal_Q(a, b, c)
```

### gslcdflognormalPinv
Link to:
```cpp
gsl_cdf_lognormal_Pinv(a, b, c)
```

### gslcdflognormalQinv
Link to:
```cpp
gsl_cdf_lognormal_Qinv(a, b, c)
```

### gslcdfgumbel1P
Link to:
```cpp
gsl_cdf_gumbel1_P(a, b, c)
```

### gslcdfgumbel1Q
Link to:
```cpp
gsl_cdf_gumbel1_Q(a, b, c)
```

### gslcdfgumbel1Pinv
Link to:
```cpp
gsl_cdf_gumbel1_Pinv(a, b, c)
```

### gslcdfgumbel1Qinv
Link to:
```cpp
gsl_cdf_gumbel1_Qinv(a, b, c)
```

### gslcdfgumbel2P
Link to:
```cpp
gsl_cdf_gumbel2_P(a, b, c)
```

### gslcdfgumbel2Q
Link to:
```cpp
gsl_cdf_gumbel2_Q(a, b, c)
```

### gslcdfgumbel2Pinv
Link to:
```cpp
gsl_cdf_gumbel2_Pinv(a, b, c)
```

### gslcdfgumbel2Qinv
Link to:
```cpp
gsl_cdf_gumbel2_Qinv(a, b, c)
```

### gslcdfweibullP
Link to:
```cpp
gsl_cdf_weibull_P(a, b, c)
```

### gslcdfweibullQ
Link to:
```cpp
gsl_cdf_weibull_Q(a, b, c)
```

### gslcdfweibullPinv
Link to:
```cpp
gsl_cdf_weibull_Pinv(a, b, c)
```

### gslcdfweibullQinv
Link to:
```cpp
gsl_cdf_weibull_Qinv(a, b, c)
```

### gslcdfparetoP
Link to:
```cpp
gsl_cdf_pareto_P(a, b, c)
```

### gslcdfparetoQ
Link to:
```cpp
gsl_cdf_pareto_Q(a, b, c)
```

### gslcdfparetoPinv
Link to:
```cpp
gsl_cdf_pareto_Pinv(a, b, c)
```

### gslcdfparetoQinv
Link to:
```cpp
gsl_cdf_pareto_Qinv(a, b, c)
```

### gslcdflogisticP
Link to:
```cpp
gsl_cdf_logistic_P(a, b)
```

### gslcdflogisticQ
Link to:
```cpp
gsl_cdf_logistic_Q(a, b)
```

### gslcdflogisticPinv
Link to:
```cpp
gsl_cdf_logistic_Pinv(a, b)
```

### gslcdflogisticQinv
Link to:
```cpp
gsl_cdf_logistic_Qinv(a, b)
```

### gslcdfbinomialP
Link to:
```cpp
gsl_cdf_binomial_P(a, b, c)
```

### gslcdfbinomialQ
Link to:
```cpp
gsl_cdf_binomial_Q(a, b, c)
```

### gslcdfpoissonP
Link to:
```cpp
gsl_cdf_poisson_P(a, b)
```

### gslcdfpoissonQ
Link to:
```cpp
gsl_cdf_poisson_Q(a, b)
```

### gslcdfgeometricP
Link to:
```cpp
gsl_cdf_geometric_P(a, b)
```

### gslcdfgeometricQ
Link to:
```cpp
gsl_cdf_geometric_Q(a, b)
```

### gslcdfnegativebinomialP
Link to:
```cpp
gsl_cdf_negative_binomial_P(a, b, c)
```

### gslcdfnegativebinomialQ
Link to:
```cpp
gsl_cdf_negative_binomial_Q(a, b, c)
```

### gslcdfpascalP
Link to:
```cpp
gsl_cdf_pascal_P(a, b, c)
```

### gslcdfpascalQ
Link to:
```cpp
gsl_cdf_pascal_Q(a, b, c)
```

### gslranbernoullipdf
Link to:
```cpp
gsl_ran_bernoulli_pdf(a, b)
```

### gslranbeta
Link to:
```cpp
gsl_ran_beta(a, b, c)
```

### gslranbetapdf
Link to:
```cpp
gsl_ran_beta_pdf(a, b, c)
```

### gslranbinomialpdf
Link to:
```cpp
gsl_ran_binomial_pdf(a, b, c)
```

### gslranexponential
Link to:
```cpp
gsl_ran_exponential(a, b)
```

### gslranexponentialpdf
Link to:
```cpp
gsl_ran_exponential_pdf(a, b)
```

### gslranexppow
Link to:
```cpp
gsl_ran_exppow(a, b, c)
```

### gslranexppowpdf
Link to:
```cpp
gsl_ran_exppow_pdf(a, b, c)
```

### gslrancauchy
Link to:
```cpp
gsl_ran_cauchy(a, b)
```

### gslrancauchypdf
Link to:
```cpp
gsl_ran_cauchy_pdf(a, b)
```

### gslranchisq
Link to:
```cpp
gsl_ran_chisq(a, b)
```

### gslranchisqpdf
Link to:
```cpp
gsl_ran_chisq_pdf(a, b)
```

### gslranerlang
Link to:
```cpp
gsl_ran_erlang(a, b, c)
```

### gslranerlangpdf
Link to:
```cpp
gsl_ran_erlang_pdf(a, b, c)
```

### gslranfdist
Link to:
```cpp
gsl_ran_fdist(a, b, c)
```

### gslranfdistpdf
Link to:
```cpp
gsl_ran_fdist_pdf(a, b, c)
```

### gslranflat
Link to:
```cpp
gsl_ran_flat(a, b, c)
```

### gslranflatpdf
Link to:
```cpp
gsl_ran_flat_pdf(a, b, c)
```

### gslrangamma
Link to:
```cpp
gsl_ran_gamma(a, b, c)
```

### gslrangammaint
Link to:
```cpp
gsl_ran_gamma_int(a, b, c)
```

### gslrangammapdf
Link to:
```cpp
gsl_ran_gamma_pdf(a, b, c)
```

### gslrangammamt
Link to:
```cpp
gsl_ran_gamma_mt(a, b, c)
```

### gslrangammaknuth
Link to:
```cpp
gsl_ran_gamma_knuth(a, b, c)
```

### gslrangaussian
Link to:
```cpp
gsl_ran_gaussian(a, b)
```

### gslrangaussianratiomethod
Link to:
```cpp
gsl_ran_gaussian_ratio_method(a, b)
```

### gslrangaussianziggurat
Link to:
```cpp
gsl_ran_gaussian_ziggurat(a, b)
```

### gslrangaussianpdf
Link to:
```cpp
gsl_ran_gaussian_pdf(a, b)
```

### gslranugaussian
Link to:
```cpp
gsl_ran_ugaussian(a)
```

### gslranugaussianratiomethod
Link to:
```cpp
gsl_ran_ugaussian_ratio_method(a)
```

### gslranugaussianpdf
Link to:
```cpp
gsl_ran_ugaussian_pdf(a)
```

### gslrangaussiantail
Link to:
```cpp
gsl_ran_gaussian_tail(a, b, c)
```

### gslrangaussiantailpdf
Link to:
```cpp
gsl_ran_gaussian_tail_pdf(a, b, c)
```

### gslranugaussiantail
Link to:
```cpp
gsl_ran_ugaussian_tail(a, b)
```

### gslranugaussiantailpdf
Link to:
```cpp
gsl_ran_ugaussian_tail_pdf(a, b)
```

### gslranlandau
Link to:
```cpp
gsl_ran_landau(a)
```

### gslranlandaupdf
Link to:
```cpp
gsl_ran_landau_pdf(a)
```

### gslrangeometricpdf
Link to:
```cpp
gsl_ran_geometric_pdf(a, b)
```

### gslrangumbel1
Link to:
```cpp
gsl_ran_gumbel1(a, b, c)
```

### gslrangumbel1pdf
Link to:
```cpp
gsl_ran_gumbel1_pdf(a, b, c)
```

### gslrangumbel2
Link to:
```cpp
gsl_ran_gumbel2(a, b, c)
```

### gslrangumbel2pdf
Link to:
```cpp
gsl_ran_gumbel2_pdf(a, b, c)
```

### gslranlogistic
Link to:
```cpp
gsl_ran_logistic(a, b)
```

### gslranlogisticpdf
Link to:
```cpp
gsl_ran_logistic_pdf(a, b)
```

### gslranlognormal
Link to:
```cpp
gsl_ran_lognormal(a, b, c)
```

### gslranlognormalpdf
Link to:
```cpp
gsl_ran_lognormal_pdf(a, b, c)
```

### gslranlogarithmicpdf
Link to:
```cpp
gsl_ran_logarithmic_pdf(a, b)
```

### gslrannegativebinomialpdf
Link to:
```cpp
gsl_ran_negative_binomial_pdf(a, b, c)
```

### gslranpascalpdf
Link to:
```cpp
gsl_ran_pascal_pdf(a, b, c)
```

### gslranpareto
Link to:
```cpp
gsl_ran_pareto(a, b, c)
```

### gslranparetopdf
Link to:
```cpp
gsl_ran_pareto_pdf(a, b, c)
```

### gslranpoissonpdf
Link to:
```cpp
gsl_ran_poisson_pdf(a, b)
```

### gslranrayleigh
Link to:
```cpp
gsl_ran_rayleigh(a, b)
```

### gslranrayleighpdf
Link to:
```cpp
gsl_ran_rayleigh_pdf(a, b)
```

### gslranrayleightail
Link to:
```cpp
gsl_ran_rayleigh_tail(a, b, c)
```

### gslranrayleightailpdf
Link to:
```cpp
gsl_ran_rayleigh_tail_pdf(a, b, c)
```

### gslrantdist
Link to:
```cpp
gsl_ran_tdsit(a, b)
```

### gslrantdistpdf
Link to:
```cpp
gsl_ran_tdsit_pdf(a, b)
```

### gslranlaplace
Link to:
```cpp
gsl_ran_laplace(a, b)
```

### gslranlaplacepdf
Link to:
```cpp
gsl_ran_laplace_pdf(a, b)
```

### gslranlevy
Link to:
```cpp
gsl_ran_levy(a, b, c)
```

### gslranweibull
Link to:
```cpp
gsl_ran_weibull(a, b, c)
```

### gslranweibullpdf
Link to:
```cpp
gsl_ran_weibull_pdf(a, b, c)
```

### gslsfairyAi
Link to:
```cpp
gsl_sf_airy_Ai(a, b)
```

### gslsfairyBi
Link to:
```cpp
gsl_sf_airy_Bi(a, b)
```

### gslsfairyAiscaled
Link to:
```cpp
gsl_sf_airy_Ai_scaled(a, b)
```

### gslsfairyBiscaled
Link to:
```cpp
gsl_sf_airy_Bi_scaled(a, b)
```

### gslsfairyAideriv
Link to:
```cpp
gsl_sf_airy_Ai_deriv(a, b)
```

### gslsfairyBideriv
Link to:
```cpp
gsl_sf_airy_Bi_deriv(a, b)
```

### gslsfairyAiderivscaled
Link to:
```cpp
gsl_sf_airy_Ai_deriv_scaled(a, b)
```

### gslsfairyBiderivscaled
Link to:
```cpp
gsl_sf_airy_Bi_deriv_scaled(a, b)
```

### gslsfairyzeroAi
Link to:
```cpp
gsl_sf_airy_Ai(a, b)
```

### gslsfairyzeroBi
Link to:
```cpp
gsl_sf_airy_aero_Bi(a)
```

### gslsfairyzeroAideriv
Link to:
```cpp
gsl_sf_airy_aero_Ai_deriv(a)
```

### gslsfairyzeroBideriv
Link to:
```cpp
gsl_sf_airy_aero_Bi_deriv(a)
```

### gslsfbesselJ0
Link to:
```cpp
gsl_sf_bessel_J0(a)
```

### gslsfbesselJ1
Link to:
```cpp
gsl_sf_bessel_J1(a)
```

### gslsfbesselJn
Link to:
```cpp
gsl_sf_bessel_Jn(a, b)
```

### gslsfbesselY0
Link to:
```cpp
gsl_sf_bessel_Y0(a)
```

### gslsfbesselY1
Link to:
```cpp
gsl_sf_bessel_Y1(a)
```

### gslsfbesselYn
Link to:
```cpp
gsl_sf_bessel_Yn(a, b)
```

### gslsfbesselI0
Link to:
```cpp
gsl_sf_bessel_I0(a)
```

### gslsfbesselI1
Link to:
```cpp
gsl_sf_bessel_I1(a)
```

### gslsfbesselIn
Link to:
```cpp
gsl_sf_bessel_In(a, b)
```

### gslsfbesselI0scaled
Link to:
```cpp
gsl_sf_bessel_I0_scaled(a)
```

### gslsfbesselI1scaled
Link to:
```cpp
gsl_sf_bessel_I1_scaled(a)
```

### gslsfbesselInscaled
Link to:
```cpp
gsl_sf_bessel_In_scaled(a, b)
```

### gslsfbesselK0
Link to:
```cpp
gsl_sf_bessel_K0(a)
```

### gslsfbesselK1
Link to:
```cpp
gsl_sf_bessel_K1(a)
```

### gslsfbesselKn
Link to:
```cpp
gsl_sf_bessel_Kn(a, b)
```

### gslsfbesselK0scaled
Link to:
```cpp
gsl_sf_bessel_K0_scaled(a)
```

### gslsfbesselK1scaled
Link to:
```cpp
gsl_sf_bessel_K1_scaled(a)
```

### gslsfbesselKnscaled
Link to:
```cpp
gsl_sf_bessel_Kn_scaled(a, b)
```

### gslsfbesselj0
Link to:
```cpp
gsl_sf_bessel_j0(a)
```

### gslsfbesselj1
Link to:
```cpp
gsl_sf_bessel_j1(a)
```

### gslsfbesselj2
Link to:
```cpp
gsl_sf_bessel_j2(a)
```

### gslsfbesseljl
Link to:
```cpp
gsl_sf_bessel_jl(a, b)
```

### gslsfbessely0
Link to:
```cpp
gsl_sf_bessel_y0(a)
```

### gslsfbessely1
Link to:
```cpp
gsl_sf_bessel_y0(a)
```

### gslsfbessely2
Link to:
```cpp
gsl_sf_bessel_y0(a)
```

### gslsfbesselyl
Link to:
```cpp
gsl_sf_bessel_jl(a, b)
```

### gslsfbesseli0scaled
Link to:
```cpp
gsl_sf_bessel_i0_scaled(a)
```

### gslsfbesseli1scaled
Link to:
```cpp
gsl_sf_bessel_i1_scaled(a)
```

### gslsfbesseli2scaled
Link to:
```cpp
gsl_sf_bessel_i2_scaled(a)
```

### gslsfbesselilscaled
Link to:
```cpp
gsl_sf_bessel_il_scaled(a, b)
```

### gslsfbesselk0scaled
Link to:
```cpp
gsl_sf_bessel_k0_scaled(a)
```

### gslsfbesselk1scaled
Link to:
```cpp
gsl_sf_bessel_k1_scaled(a)
```

### gslsfbesselk2scaled
Link to:
```cpp
gsl_sf_bessel_k2_scaled(a)
```

### gslsfbesselklscaled
Link to:
```cpp
gsl_sf_bessel_kl_scaled(a, b)
```

### gslsfbesselJnu
Link to:
```cpp
gsl_sf_bessel_Jnu(a, b)
```

### gslsfbesselYnu
Link to:
```cpp
gsl_sf_bessel_Ynu(a, b)
```

### gslsfbesselInuscaled
Link to:
```cpp
gsl_sf_bessel_Inu_scaled(a, b)
```

### gslsfbesselInu
Link to:
```cpp
gsl_sf_bessel_Inu(a, b)
```

### gslsfbesselKnuscaled
Link to:
```cpp
gsl_sf_bessel_Knu_scaled(a, b)
```

### gslsfbesselKnu
Link to:
```cpp
gsl_sf_bessel_Knu(a, b)
```

### gslsfbessellnKnu
Link to:
```cpp
gsl_sf_bessel_lnKnu(a, b)
```

### gslsfbesselzeroJ0
Link to:
```cpp
gsl_sf_bessel_zero_J0(a)
```

### gslsfbesselzeroJ1
Link to:
```cpp
gsl_sf_bessel_zero_J1(a)
```

### gslsfbesselzeroJnu
Link to:
```cpp
gsl_sf_bessel_zero_Jnu(a, b)
```

### gslsfclausen
Link to:
```cpp
gsl_sf_clausen(a)
```

### gslsfhydrogenicR1
Link to:
```cpp
gsl_sf_hydrogenicR_1(a, b)
```

### gslsfdawson
Link to:
```cpp
gsl_sf_dawnson(a)
```

### gslsfdebye1
Link to:
```cpp
gsl_sf_debye_1(a)
```

### gslsfdebye2
Link to:
```cpp
gsl_sf_debye_2(a)
```

### gslsfdebye3
Link to:
```cpp
gsl_sf_debye_3(a)
```

### gslsfdebye4
Link to:
```cpp
gsl_sf_debye_4(a)
```

### gslsfdebye5
Link to:
```cpp
gsl_sf_debye_5(a)
```

### gslsfdebye6
Link to:
```cpp
gsl_sf_debye_6(a)
```

### gslsfdilog
Link to:
```cpp
gsl_sf_dilog(a)
```

### gslsfmultiply
Link to:
```cpp
gsl_sf_multiply(a, b)
```

### gslsfellintKcomp
Link to:
```cpp
gsl_sf_ellint_Kcomp(a, b)
```

### gslsfellintEcomp
Link to:
```cpp
gsl_sf_ellint_Ecomp(a, b)
```

### gslsfellintPcomp
Link to:
```cpp
gsl_sf_ellint_Pcomp(a, b, c)
```

### gslsfellintDcomp
Link to:
```cpp
gsl_sf_ellint_Dcomp(a, b)
```

### gslsfellintF
Link to:
```cpp
gsl_sf_ellint_F(a, b, c)
```

### gslsfellintE
Link to:
```cpp
gsl_sf_ellint_E(a, b, c)
```

### gslsfellintRC
Link to:
```cpp
gsl_sf_ellint_RC(a, b, c)
```

### gslsferfc
Link to:
```cpp
gsl_sf_erfc(a)
```

### gslsflogerfc
Link to:
```cpp
gsl_sf_log_erfc(a)
```

### gslsferf
Link to:
```cpp
gsl_sf_erf(a)
```

### gslsferfZ
Link to:
```cpp
gsl_sf_erf_Z(a)
```

### gslsferfQ
Link to:
```cpp
gsl_sf_erf_Q(a)
```

### gslsfhazard
Link to:
```cpp
gsl_sf_hazard(a)
```

### gslsfexp
Link to:
```cpp
gsl_sf_exp(a)
```

### gslsfexpmult
Link to:
```cpp
gsl_sf_exp_mult(a, b)
```

### gslsfexpm1
Link to:
```cpp
gsl_sf_expm1(a)
```

### gslsfexprel
Link to:
```cpp
gsl_sf_exprel(a)
```

### gslsfexprel2
Link to:
```cpp
gsl_sf_exprel_2(a)
```

### gslsfexpreln
Link to:
```cpp
gsl_sf_exprel_n(a, b)
```

### gslsfexpintE1
Link to:
```cpp
gsl_sf_expint_E1(a)
```

### gslsfexpintE2
Link to:
```cpp
gsl_sf_expint_E2(a)
```

### gslsfexpintEn
Link to:
```cpp
gsl_sf_expint_En(a, b)
```

### gslsfexpintE1scaled
Link to:
```cpp
gsl_sf_expint_E1_scaled(a)
```

### gslsfexpintE2scaled
Link to:
```cpp
gsl_sf_expint_E1_scaled(a)
```

### gslsfexpintEnscaled
Link to:
```cpp
gsl_sf_expint_En_scaled(a, b)
```

### gslsfexpintEi
Link to:
```cpp
gsl_sf_expint_Ei(a)
```

### gslsfexpintEiscaled
Link to:
```cpp
gsl_sf_expint_Ei_scaled(a)
```

### gslsfShi
Link to:
```cpp
gsl_sf_Shi(a)
```

### gslsfChi
Link to:
```cpp
gsl_sf_Chi(a)
```

### gslsfexpint3
Link to:
```cpp
gsl_sf_expint_3(a)
```

### gslsfSi
Link to:
```cpp
gsl_sf_Si(a)
```

### gslsfCi
Link to:
```cpp
gsl_sf_Ci(a)
```

### gslsfatanint
Link to:
```cpp
gsl_sf_atanint(a)
```

### gslsffermidiracm1
Link to:
```cpp
gsl_sf_fermi_dirac_m1(a)
```

### gslsffermidirac0
Link to:
```cpp
gsl_sf_fermi_dirac_0(a)
```

### gslsffermidirac1
Link to:
```cpp
gsl_sf_fermi_dirac_1(a)
```

### gslsffermidirac2
Link to:
```cpp
gsl_sf_fermi_dirac_2(a)
```

### gslsffermidiracint
Link to:
```cpp
gsl_sf_fermi_dirac_int(a, b)
```

### gslsffermidiracmhalf
Link to:
```cpp
gsl_sf_fermi_dirac_mhalf(a)
```

### gslsffermidirachalf
Link to:
```cpp
gsl_sf_fermi_dirac_half(a)
```

### gslsffermidirac3half
Link to:
```cpp
gsl_sf_fermi_dirac_3half(a)
```

### gslsffermidiracinc0
Link to:
```cpp
gsl_sf_fermi_dirac_inc_0(a, b)
```

### gslsflngamma
Link to:
```cpp
gsl_sf_lngamma(a)
```

### gslsfgamma
Link to:
```cpp
gsl_sf_gamma(a)
```

### gslsfgammastar
Link to:
```cpp
gsl_sf_gammastar(a)
```

### gslsfgammainv
Link to:
```cpp
gsl_sf_gammainv(a)
```

### gslsftaylorcoeff
Link to:
```cpp
gsl_sf_taylorcoeff(a, b)
```

### gslsffact
Link to:
```cpp
gsl_sf_fact(a)
```

### gslsfdoublefact
Link to:
```cpp
gsl_sf_doublefact(a)
```

### gslsflnfact
Link to:
```cpp
gsl_sf_lnfact(a)
```

### gslsflndoublefact
Link to:
```cpp
gsl_sf_lndoublefact(a)
```

### gslsflnchoose
Link to:
```cpp
gsl_sf_lnchoose(a, b)
```

### gslsfchoose
Link to:
```cpp
gsl_sf_choose(a, b)
```

### gslsflnpoch
Link to:
```cpp
gsl_sf_lnpoch(a, b)
```

### gslsfpoch
Link to:
```cpp
gsl_sf_poch(a, b)
```

### gslsfpochrel
Link to:
```cpp
gsl_sf_pochrel(a, b)
```

### gslsfgammaincQ
Link to:
```cpp
gsl_sf_gamma_inc_Q(a, b)
```

### gslsfgammaincP
Link to:
```cpp
gsl_sf_gamma_inc_P(a, b)
```

### gslsfgammainc
Link to:
```cpp
gsl_sf_gamma_inc(a, b)
```

### gslsflnbeta
Link to:
```cpp
gsl_sf_lnbeta(a, b)
```

### gslsfbeta
Link to:
```cpp
gsl_sf_beta(a, b)
```

### gslsfbetainc
Link to:
```cpp
gsl_sf_betaçinc(a, b, c)
```

### gslsfgegenpoly1
Link to:
```cpp
gsl_sf_gegenpoly_1(a, b)
```

### gslsfgegenpoly2
Link to:
```cpp
gsl_sf_gegenpoly_2(a, b)
```

### gslsfgegenpoly3
Link to:
```cpp
gsl_sf_gegenpoly_3(a, b)
```

### gslsfgegenpolyn
Link to:
```cpp
gsl_sf_gegenpoly_n(a, b, c)
```

### gslsfhyperg0F1
Link to:
```cpp
gsl_sf_hyperg_0F1(a, b)
```

### gslsfhyperg1F1int
Link to:
```cpp
gsl_sf_hyperg_1F1_inc(a, b, c)
```

### gslsfhyperg1F1
Link to:
```cpp
gsl_sf_hyperg_1F1(a, b, c)
```

### gslsfhypergUint
Link to:
```cpp
gsl_sf_hyperg_U_inc(a, b, c)
```

### gslsfhypergU
Link to:
```cpp
gsl_sf_hyperg_U(a, b, c)
```

### gslsfhyperg2F0
Link to:
```cpp
gsl_sf_hyperg_U_2F0(a, b, c)
```

### gslsflaguerre1
Link to:
```cpp
gsl_sf_laguerre_1(a, b)
```

### gslsflaguerre2
Link to:
```cpp
gsl_sf_laguerre_2(a, b)
```

### gslsflaguerre3
Link to:
```cpp
gsl_sf_laguerre_3(a, b)
```

### gslsflaguerren
Link to:
```cpp
gsl_sf_laguerre_n(a, b, c)
```

### gslsflambertW0
Link to:
```cpp
gsl_sf_lambert_W0(a)
```

### gslsflambertWm1
Link to:
```cpp
gsl_sf_lambert_Wm1(a)
```

### gslsflegendrePl
Link to:
```cpp
gsl_sf_legendre_Pl(a, b)
```

### gslsflegendreP1
Link to:
```cpp
gsl_sf_legendre_P1(a)
```

### gslsflegendreP2
Link to:
```cpp
gsl_sf_legendre_P2(a)
```

### gslsflegendreP3
Link to:
```cpp
gsl_sf_legendre_P3(a)
```

### gslsflegendreQ0
Link to:
```cpp
gsl_sf_legendre_Q0(a)
```

### gslsflegendreQ1
Link to:
```cpp
gsl_sf_legendre_Q1(a)
```

### gslsflegendreQl
Link to:
```cpp
gsl_sf_legendre_Ql(a, b)
```

### gslsflegendrePlm
Link to:
```cpp
gsl_sf_legendre_Plm(a, b, c)
```

### gslsflegendresphPlm
Link to:
```cpp
gsl_sf_legendre_sphP1m(a, b, c)
```

### gslsflegendrearraysize
Link to:
```cpp
gsl_sf_legendre_array_size(a, b)
```

### gslsfconicalPhalf
Link to:
```cpp
gsl_sf_conicalP_half(a, b)
```

### gslsfconicalPmhalf
Link to:
```cpp
gsl_sf_conicalP_mhalf(a, b)
```

### gslsfconicalP0
Link to:
```cpp
gsl_sf_conicalP_0(a, b)
```

### gslsfconicalP1
Link to:
```cpp
gsl_sf_conicalP_1(a, b)
```

### gslsfconicalPsphreg
Link to:
```cpp
gsl_sf_conicalP_sph_reg(a, b, c)
```

### gslsfconicalPcylreg
Link to:
```cpp
gsl_sf_conicalP_cyl_reg(a, b, c)
```

### gslsflegendreH3d0
Link to:
```cpp
gsl_sf_legendre_H3d_0(a, b)
```

### gslsflegendreH3d1
Link to:
```cpp
gsl_sf_legendre_H3d_1(a, b)
```

### gslsflegendreH3d
Link to:
```cpp
gsl_sf_legendre_H3d(a, b, c)
```

### gslsflog
Link to:
```cpp
gsl_sf_log(a)
```

### gslsflogabs
Link to:
```cpp
gsl_sf_log_abs(a)
```

### gslsflog1plusx
Link to:
```cpp
gsl_sf_log_1plusx(a)
```

### gslsflog1plusxmx
Link to:
```cpp
gsl_sf_log_1plusx_mx(a)
```

### gslsfpowint
Link to:
```cpp
gsl_sf_pow_int(a, b)
```

### gslsfpsiint
Link to:
```cpp
gsl_sf_psi_int(a)
```

### gslsfpsi
Link to:
```cpp
gsl_sf_psi(a)
```

### gslsfpsi1piy
Link to:
```cpp
gsl_sf_psi_1piy(a)
```

### gslsfpsi1int
Link to:
```cpp
gsl_sf_psi_1_int(a)
```

### gslsfpsi1
Link to:
```cpp
gsl_sf_psi_1(a)
```

### gslsfpsin
Link to:
```cpp
gsl_sf_psi_n(a, b)
```

### gslsfsynchrotron1
Link to:
```cpp
gsl_sf_synchrotron_1(a)
```

### gslsfsynchrotron2
Link to:
```cpp
gsl_sf_synchrotron_2(a)
```

### gslsftransport2
Link to:
```cpp
gsl_sf_transport_2(a)
```

### gslsftransport3
Link to:
```cpp
gsl_sf_transport_3(a)
```

### gslsftransport4
Link to:
```cpp
gsl_sf_transport_4(a)
```

### gslsftransport5
Link to:
```cpp
gsl_sf_transport_5(a)
```

### gslsfsin
Link to:
```cpp
gsl_sf_sin(a)
```

### gslsfcos
Link to:
```cpp
gsl_sf_cos(a)
```

### gslsfhypot
Link to:
```cpp
gsl_sf_hypot(a, b)
```

### gslsfsinc
Link to:
```cpp
gsl_sf_sinc(a)
```

### gslsflnsinh
Link to:
```cpp
gsl_sf_lnsinh(a)
```

### gslsflncosh
Link to:
```cpp
gsl_sf_lncosh(a)
```

### gslsfanglerestrictsymm
Link to:
```cpp
gsl_sf_andle_restrict_symm(a)
```

### gslsfanglerestrictpos
Link to:
```cpp
gsl_sf_angle_restrict_pos(a)
```

### gslsfzetaint
Link to:
```cpp
gsl_sf_zeta_int(a)
```

### gslsfzeta
Link to:
```cpp
gsl_sf_zeta(a)
```

### gslsfzetam1
Link to:
```cpp
gsl_sf_zetam1(a)
```

### gslsfzetam1int
Link to:
```cpp
gsl_sf_zetam1_int(a)
```

### gslsfhzeta
Link to:
```cpp
gsl_sf_hzeta(a, b)
```

### gslsfetaint
Link to:
```cpp
gsl_sf_eta_int(a)
```

### gslsfeta
Link to:
```cpp
gsl_sf_eta(a)
```

## ff-Ipopt

Refer to the [Ipopt documentation](https://projects.coin-or.org/Ipopt) for more informations.

### IPOPT

$\codered$

## fflapack

Refer to the [LAPACK documentation](http://www.netlib.org/lapack/) for more informations.

###inv

$\codered$

### dgeev

$\codered$

### zgeev

$\codered$

### geev

$\codered$

### geev

$\codered$

### dggev

$\codered$

### zggev

$\codered$

### dsygvd

$\codered$

### dgesdd

$\codered$

### zhegv

$\codered$

### dsyev

$\codered$

### zheev

$\codered$

## ff-mmap-semaphore

### Wait

$\codered$

### trywait

$\codered$

### Post

$\codered$

### msync

$\codered$

### Read

$\codered$

### Write

$\codered$

## ffnewuoa

### newuoa

$\codered$

## ff-NLopt

Refer to the [NLOPT documentation](https://nlopt.readthedocs.io/en/latest/) for more informations.

### nloptDIRECT

$\codered$

### nloptDIRECTL

$\codered$

### nloptDIRECTLRand

$\codered$

### nloptDIRECTScal

$\codered$

### nloptDIRECTNoScal

$\codered$

### nloptDIRECTLNoScal

$\codered$

### nloptDIRECTLRandNoScal

$\codered$

### nloptOrigDIRECT

$\codered$

### nloptOrigDIRECTL

$\codered$

### nloptStoGO

$\codered$

### nloptStoGORand

$\codered$

### nloptLBFGS

$\codered$

### nloptPRAXIS

$\codered$

### nloptVar1

$\codered$

### nloptVar2

$\codered$

### nloptTNewton

$\codered$

### nloptTNewtonRestart

$\codered$

### nloptTNewtonPrecond

$\codered$

### nloptTNewtonPrecondRestart

$\codered$

### nloptCRS2

$\codered$

### nloptMMA

$\codered$

### nloptCOBYLA

$\codered$

### nloptNEWUOA

$\codered$

### nloptNEWUOABound

$\codered$

### nloptNelderMead

$\codered$

### nloptSbplx

$\codered$

### nloptBOBYQA

$\codered$

### nloptISRES

$\codered$

### nloptSLSQP

$\codered$

### nloptMLSL

$\codered$

### nloptMLSLLDS

$\codered$

### nloptAUGLAG

$\codered$

### nloptAUGLAGEQ

$\codered$

## ffrandom

### srandomdev

$\codered$

### srandom

$\codered$

### random

$\codered$

## FreeFemQA

### MeshGenQA

$\codered$

## freeyams

### freeyams

$\codered$

<!---
## funcTempate

$\codered$
--->

## gmsh

Need
```freefem
load "gsmh"
```

The `gmsh` software is available [here](http://gmsh.info/)

### gmshload

Load a 2D mesh build with Gmsh.

```freefem
mesh Th = gmshload(MeshFile, [reftri=RefTri], [renum=Renum]);
```

<u>Parameters:</u>

- `MeshFile` (`:::freefem string`)<br/>
Mesh file name
- _`:::freefem reftri=`_ ($\codered$)
- _`:::freefem renum=`_ ($\codered$)

<u>Output:</u>

- `Th` (`:::freefem mesh`)

### gmshload3

Load a 3D mesh build with Gmsh.

```freefem
mesh3 Th = gmshload(MeshFile, [reftet=RefTet], [renum=Renum]);
```

<u>Parameters:</u>

- `MeshFile` (`:::freefem string`)<br/>
Mesh file name
- _`:::freefem reftet=`_ ($\codered$)
- _`:::freefem renum=`_ ($\codered$)

<u>Output:</u>

- `Th` (`:::freefem mesh3`)

### savegmsh

$\codered$

## gsl

### gslpolysolvequadratic

$\codered$

### gslpolysolvecubic

$\codered$

### gslpolycomplexsolve

$\codered$

### gslrnguniform

$\codered$

### gslrnguniformpos

$\codered$

### gslname

$\codered$

### gslrngget

$\codered$

### gslrngmin

$\codered$

### gslrngmax

$\codered$

### gslrngset

$\codered$

### gslrngtype

$\codered$

## ilut

### applyIlutPrecond

$\codered$

### makeIlutPrecond

$\codered$

## iohdf5

### savehdf5sol

$\codered$

## iovtk

### savevtk

Save mesh or solution in vtk/vtu format.

```freefem
savetk(FileName, Th, [Ux, Uy, Uz], p, [dataname=DataName], [withsurfacemesh=WithSurfaceMesh], [order=Order], [floatmesh=FloatMesh], [floatsol=FloatSol], [bin=Bin], [swap=Swap]);
```

<u>Parameters</u>:

 - `FileName` (`:::freefem string`)<br/>
 File name: `*.vtk` or `*.vtu`
 - `Th` (`:::freefem mesh` or `:::freefem mesh3`)
 - `[Ux, Uy, Uz], p` (`:::freefem fespace` function of vector of `:::freefem fespace` functions)<br/>
 Solutions to save, as much as wanted
 - `:::frefem dataname=` (`:::freefem string`)<br/>
 Name of solutions, seprated by a space
 - `:::freefem withsurfacemesh=` (`:::freefem bool`)<br/>
 $\codered$
 - `:::freefem order=` (`:::freefem int[int]`)</br>
 Order of solutions. Available: 0 or 1
 - `:::freefem floatmesh=` (`:::freefem bool`)<br/>
 $\codered$
 - `:::freefem floatsol=` (`:::freefem bool`)<br/>
 $\codered$
 - `:::freefem bin=` (`:::freefem bool`)<br/>
 If true, save file in binary format
 - `:::freefem swap` (`:::freefem bool`)</br>
 $\codered$

<u>Output</u>:

 - None

### vtkload

$\codered$

### vtkload3

$\codered$

## isoline

Need
```freefem
load "isoline"
```

### isoline

```freefem
int N = isoline(Th, u, xy, iso=Iso, close=Close, smoothing=Smoothing, ratio=Ratio, eps=Eps, beginend=BeginEnd, file=File);
```

$\codered$

### Curve

$\codered$

### Area

$\codered$

### findallocalmin

$\codered$

## lapack

### inv

$\codered$

### dgeev

$\codered$

### zgeev

$\codered$

### geev

$\codered$

### dggev

$\codered$

### zggev

$\codered$

### dsygvd

$\codered$

### dgesdd

$\codered$

### zhegv

$\codered$

### dsyev

$\codered$

### zheev

$\codered$

### dgelsy

$\codered$

## lgbmo

### bmo

$\codered$

## mat_dervieux

### MatUpWind1

$\codered$

## mat_psi

### MatUpWind0

$\codered$

## medit

### medit

$\codered$

### savesol

$\codered$

### readsol

$\codered$

## metis

### metisnodal

$\codered$

### metisdual

$\codered$

## MetricKuate

### MetricKuate

$\codered$

## MetricPk

### MetricPk

$\codered$

## mmg3d

### mmg3d

$\codered$

## mmg3d-v4.0

### mmg3d

$\codered$

## msh3

### change

$\codered$

### movemesh23

$\codered$

### movemesh2D3Dsurf

$\codered$

### movemesh3

$\codered$

### movemesh

$\codered$

### movemesh3D

$\codered$

### deplacement

$\codered$

### checkbemesh

$\codered$

### buildlayers

$\codered$

### bcube

$\codered$

### cube
Construct a cubic mesh.

```freefem
mesh3 Th = cube(nnX, nnY, nnZ, [X(x), Y(y), Z(z)], [label=Label], [flags=Flags], [region=Region]);
```

<u>Parameters:</u>

 - `nnX` (`:::freefem int`)<br/>
 Number of discretization point along $x$
 - `nnY` (`:::freefem int`)<br/>
 Number of discretization point along $y$
 - `nnZ` (`:::freefem int`)<br/>
  Number of discretization point along $z$
  - `X(x)` (`:::freefem func`) _[Optional]_<br/>
  Affine function of $x$ to define the length<br/>
  Default: `x`
  - `Y(y)` (`:::freefem func`) _[Optional]_<br/>
  Affine function of $y$ to define the width<br/>
  Default: `y`
  - `Z(z)` (`:::freefem func`) _[Optional]_<br/>
  Affine function of $z$ to define the height<br/>
  Default: `z`
  - _`:::freefem label=`_ (`:::freefem int[int]`) _[Optional]_<br/>
  List of surface labels<br/>
  Default: `[1, 2, 3, 4, 5, 6]`
  - _`:::freefem flags=`_ (`:::freefem int`) _[Optional]_<br/>
  Refer to [square](#square)
  - _`:::freefem region=`_ (`:::freefem int`) _[Optional]_<br/>
  Region number of the cube volume
  Default: `0`

<u>Output:</u>

 - `Th` (`:::freefem mesh3`)<br/>
 Cube mesh

### trunc

$\codered$

### gluemesh

$\codered$

### extract

$\codered$

### showborder

$\codered$

### getborder

$\codered$

### AddLayers

$\codered$

## mshmet

### mshmet

$\codered$

## MUMPS

### defaulttoMUMPSseq

$\codered$

## MUMPS_seq

### defaulttoMUMPSseq

$\codered$

<!---
## myfunction2

$\codered$
--->

## netgen

### netg

$\codered$

### netgstl

$\codered$

### netgload

$\codered$

## NewSolver

### defaulttoUMFPACK

$\codered$

## PARDISO

### defaulttoPARDISO

$\codered$

### ompsetnumthreads

$\codered$

### ompgetnumthreads

$\codered$

### ompgetmaxthreads

$\codered$

## pcm2rnm

### readpcm

$\codered$

## pipe

### flush

$\codered$

### sleep

$\codered$

### usleep

$\codered$

## qf11to25

### QF1d

$\codered$

### QF2d

$\codered$

### QF3d

$\codered$

### tripleQF

## scotch

### scotch

$\codered$

## shell

### readdir

$\codered$

### unlink

$\codered$

### rmdir

$\codered$

### cddir

$\codered$

### chdir

$\codered$

### basename

$\codered$

### dirname

$\codered$

### mkdir

$\codered$

### chmod

$\codered$

### cpfile

$\codered$

### stat

$\codered$

### isdir

$\codered$

### getenv

$\codered$

### setenv

$\codered$

### unsetenv

$\codered$

## splitedges

### SplitedgeMesh

$\codered$

## splitmesh12

### splitmesh12

$\codered$

## splitmesh3

### splitmesh3

$\codered$

## splitmesh4

### splimesh4

$\codered$

## splitmesh6

### splitmesh6

$\codered$

## SuperLu

### defaulttoSuperLu

$\codered$

## symmetrizeCSR

### symmetrizeCSR

$\codered$

## tetgen

Refer to the [Tetgen documentation](http://wias-berlin.de/software/tetgen/) for more informations.

### tetgconvexhull

$\codered$

### tetgtransfo

$\codered$

### tetg
Build a 3D mesh from a surface.

```freefem
mesh3 Th = tetg(Th0, [reftet=RefTet], [label=Label], [switch=Switch], [nbofholes=NbOfHoles], [holelist=HoleList], [nbofregions=NbOfRegions], [regionlist=RegionList], [nboffacetcl=NbOfFaceTcl], [facetcl=FaceTcl])
```

$\codered$

### tetgreconstruction

$\codered$

## UMFPACK64

### defaulttoUMFPACK64

$\codered$

## VTK_writer_3d

### Vtkaddmesh

$\codered$

###Vtkaddscalar

$\codered$

## VTK_writer

### Vtkaddmesh

$\codered$

### Vtkaddscalar
