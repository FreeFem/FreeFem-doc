#!/bin/bash

src='../docs/'
target='./docs/'
introductionPages=(
	index
	introduction/download
	introduction/contributing
	introduction/citation
	introduction/authors
)
documentationPages=(
	documentation/index
	documentation/Notations
	documentation/MeshGeneration
	documentation/FiniteElement
	documentation/Visualization
	documentation/AlgorithmsOptimization
	documentation/Parallelization
	documentation/Plugins
	documentation/Developers
)
tutorialsPages=(
	tutorials/index
	tutorials/Poisson
	tutorials/EquationsClassification
	tutorials/Membrane
	tutorials/Acoustics
	tutorials/ThermalConduction
	tutorials/FanBlade
	tutorials/RotatingHill
	tutorials/Elasticity
	tutorials/Stokes
	tutorials/NavierStokesProjection
	tutorials/NavierStokesNewton
	tutorials/ALargeFluidProblem
	tutorials/ComplexNumbers
	tutorials/OptimalControl
	tutorials/FlowWithShocks
	tutorials/HeatEquationOptimization
	tutorials/TimeDependentStokes
	tutorials/WifiPropagation
	tutorials/MatlabOctavePlot
)

## Copy src
rm -R $target
cp -R $src $target

## Convert Markdown to LaTeX function
function md2latex() {
	array=("$@")
  for md in "${array[@]}"
  do
    echo "convert $target$md"
    pandoc $target$md.md -o $target$md.tex --filter pandoc-latex-admonition
  done
}

function input2latex() {
  array=("$@")
  for tex in "${array[@]}"
  do
    echo "\include{$target$tex}" >> FreeFem++-doc.tex
  done
}

## Start conversion
md2latex "${introductionPages[@]}"
md2latex "${documentationPages[@]}"
md2latex "${tutorialsPages[@]}"

## Create entire documentation
cat header.tex > FreeFem++-doc.tex
echo "\part{Introduction}" >> FreeFem++-doc.tex
input2latex "${introductionPages[@]}"
echo "\part{Documentation}" >> FreeFem++-doc.tex
input2latex "${documentationPages[@]}"
echo "\part{Tutorials}" >> FreeFem++-doc.tex
input2latex "${tutorialsPages[@]}"
cat footer.tex >> FreeFem++-doc.tex

## LaTeX -> PDF
pdflatex -interaction=nonstopmode FreeFem++-doc.tex
pdflatex -interaction=nonstopmode FreeFem++-doc.tex
