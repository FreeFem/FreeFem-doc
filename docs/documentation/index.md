# Documentation

## Finite elements
(to merge with chap 6)

Currently FreeFem++ implements the following elements in 2d, (see section 6\ref{finite elements} for the full description)

`P0` piecewise constant,  
`P1` continuous piecewise linear,  
`P2` continuous piecewise quadratic,  
`P3` continuous piecewise cubic (need  `load "Element_P3"`).  
`P4` continuous piecewise quartic (need  `load "Element_P4"`).  
`RT0` Raviart-Thomas piecewise constant,  
`RT1` Raviart-Thomas degree 1  piecewise constant (need  `load "Element_Mixte"`).  
`BDM1` Brezzi-Douglas-Marini degree 1  piecewise constant (need  `load "Element_Mixte"`).     
`RT0Ortho` Nedelec type 1 degree 0  piecewise constant.  
`RT1Ortho` Nedelec type 1 degree 1  piecewise constant (need `load "Element_Mixte"`).  
`BDM1Ortho` Brezzi-Douglas-Marini degree 1  piecewise constant (need `load "Element_Mixte"`).  
`P1nc` piecewise linear non-conforming,  
`P1dc` piecewise linear discontinuous,  
`P2dc` piecewise quadratic discontinuous,  
`P2h` quadratic homogene continuous (without `P1`).  
`P3dc` piecewise cubic discontinuous (need `load "Element_P3dc"`).  
`P4dc` piecewise quartic discontinuous (need `load "Element_P4dc"`).  
`P1b` piecewise linear continuous plus bubble,  
`P2b` piecewise quadratic continuous plus bubble.  
`Morley` Morley finite element (need `load "Morley"`).  
`HCT` Hsieh-Clough-Tocher $C^1$ finite element (need `load "Element_HCT"` version 3.40).  
`P2BR` P2 Bernardi-Raugel finite element (need `load "BernadiRaugel.cpp"`).  
`P0edge` a finite element constant per edge  
`P1edge` to `P5edge`  a finite element polynomial  on edge (need `load "Element_PkEdge"`)  
...

Currently FreeFem++ implements the following elements in 3d, (see section \ref{finite elements} for the full description)

`P03d` piecewise constant,
`P13d` continuous piecewise linear,
`P23d` continuous piecewise quadratic,
`RT03d` Raviart-Thomas piecewise constant,
`Edge03d`,`Edge13d`,`Edge23d` The Nedelec Edge element 0,1,2 
`P1b3d` piecewise linear continuous plus bubble,
...

To get the full list, in a unix terminal, in directory **examples++-tutorial** do

```
FreeFem++ dumptable.edp
grep TypeOfFE lestables
```

Note that other elements can be added fairly easily.
