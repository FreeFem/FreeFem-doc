$$
\newcommand{\boldx}{\mathbf{x}}
\newcommand{\boldxi}{\boldsymbol{\xi}}
$$

The quadrature formula is like the following:
$$
\int_{D}{f(\boldx)} \approx \sum_{\ell=1}^{L}{\omega_\ell f(\boldxi_\ell)}
$$

## int1d

Quadrature formula on an edge.

### Notations

$|D|$ is the measure of the edge $D$.

For a shake of simplicity, we denote:
$$
f(\boldx) = g(t)
$$
with $0\leq t\leq 1$; $\boldx=(1-t)\boldx_0+t\boldx_1$.

### qf1pE
```freefem
int1d(Th, qfe=qf1pE)( ... )
```
or
```freefem
int1d(Th, qfe=qforder=2)( ... )
```

This quadrature formula is exact on $\P_1$.

$$
\int_{D}{f(\boldx)} \approx |D|g\left(\frac{1}{2}\right)
$$

### qf2pE
```freefem
int1d(Th, qfe=qf2pE)( ... )
```
or
```freefem
int1d(Th, qfe=qforder=3)( ... )
```

This quadrature formula is exact on $\P_3$.

$$
\int_{D}{f(\boldx)} \approx \frac{|D|}{2}\left(
	  g\left( \frac{1+\sqrt{1/3}}{2} \right)
	+ g\left( \frac{1-\sqrt{1/3}}{2} \right)
\right)
$$

### qf3pE
```freefem
int1d(Th, qfe=qf3pE)( ... )
```
or
```freefem
int1d(Th, qfe=qforder=6)( ... )
```

This quadrature formula is exact on $\P_5$.

$$
\int_{D}{f(\boldx)} \approx \frac{|D|}{18}\left(
	  5g\left( \frac{1+\sqrt{3/5}}{2} \right)
	+ 8g\left( \frac{1}{2} \right)
	+ 5g\left( \frac{1-\sqrt{3/5}}{2} \right)
\right)
$$

### qf4pE
```freefem
int1d(Th, qfe=qf4pE)( ... )
```
or
```freefem
int1d(Th, qfe=qforder=8)( ... )
```

This quadrature formula is exact on $\P_7$.

$$
\int_{D}{f(\boldx)} \approx \frac{|D|}{72}\left(
	  (18-\sqrt{30})g\left( \frac{1-\frac{\sqrt{525+70\sqrt{30}}}{35}}{2} \right)
	+ (18-\sqrt{30})g\left( \frac{1+\frac{\sqrt{525+70\sqrt{30}}}{35}}{2} \right)
	+ (18+\sqrt{30})g\left( \frac{1-\frac{\sqrt{525-70\sqrt{30}}}{35}}{2} \right)
	+ (18+\sqrt{30})g\left( \frac{1+\frac{\sqrt{525-70\sqrt{30}}}{35}}{2} \right)
\right)
$$

### qf5pE
```freefem
int1d(Th, qfe=qf5pE)( ... )
```
or
```freefem
int1d(Th, qfe=qforder=10)( ... )
```

This quadrature formula is exact on $\P_9$.

$$
\int_{D}{f(\boldx)} \approx |D|\left(
	  \frac{(332-13\sqrt{70})}{1800}g\left( \frac{1-\frac{\sqrt{245+14\sqrt{70}}}{21}}{2} \right)
	+ \frac{(332-13\sqrt{70})}{1800}g\left( \frac{1+\frac{\sqrt{245+14\sqrt{70}}}{21}}{2} \right)
	+ \frac{64}{225}g\left( \frac{1}{2} \right)
	+ \frac{(332+13\sqrt{70})}{1800}g\left( \frac{1-\frac{\sqrt{245-14\sqrt{70}}}{21}}{2} \right)
	+ \frac{(332+13\sqrt{70})}{1800}g\left( \frac{1+\frac{\sqrt{245-14\sqrt{70}}}{21}}{2} \right)
\right)
$$


### qf1pElump
```freefem
int1d(Th, qfe=qf1pElump)( ... )
```

This quadrature formula is exact on $\P_2$.

$$
\int_{D}{f(\boldx)} \approx \frac{|D|}{2}\left(
	  g\left( 0 \right)
	+ g\left( 1 \right)
\right)
$$

## int2d

### qf1pT

$\codered$

### qf2pT

$\codered$

### qf5pT

$\codered$

### qf1pTlump

$\codered$

### qf2pT4P1

$\codered$

### qf7pT

$\codered$

### qf9pT

$\codered$

## int3d

### qfV1

$\codered$

### qfV2

$\codered$

### qfV5

$\codered$

### qfV1lump

$\codered$
