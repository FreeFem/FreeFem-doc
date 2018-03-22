## documentation directory
docDirectory = "../docs/"

## list of files by category
#index
indexFilesList = [
	"index.md"
	]
indexNamesList = [
	"Home"
	]
indexTODOFile = "TODO.md"

#documentation
documentationFilesList = [
	"documentation/index.md",
	"documentation/MeshGeneration.md",
	"documentation/FiniteElement.md",
	"documentation/Visualization.md"
	]
documentationNamesList = [
	"Home",
	"Mesh generation",
	"Finite element",
	"Visualization"
	]
documentationTODOFile = "documentation/TODO.md"

#examples
examplesFilesList = [
	"examples/index.md"
	]
examplesNamesList = [
	"Home"
	]
examplesTODOFile = "examples/TODO.md"

#introduction
introductionFilesList = [
	"introduction/download.md",
	"introduction/installation.md",
	"introduction/Contributing.md"
	]
introductionNamesList = [
	"Download",
	"Installation",
	"How to contribute?"
	]
introductionTODOFile = "introduction/TODO.md"

#reference
referenceFilesList = [
	"reference/index.md",
	"reference/Types.md",
	"reference/GlobalVariables.md",
	"reference/QuadratureFormulae.md",
	"reference/Operators.md",
	"reference/Loops.md",
	"reference/IO.md",
	"reference/Functions.md",
	"reference/ExternalLibraries.md"
	]
#"reference/examples.md"
referenceNamesList = [
	"Home",
	"Types",
	"Global variables",
	"Quadrature formulae",
	"Operators",
	"Loops",
	"I/O",
	"Functions",
	"External libraries"
]
#"Examples"
referenceTODOFile = "reference/TODO.md"

#tutorial
tutorialFilesList = [
	"tutorial/index.md",
	"tutorial/Poisson.md",
	"tutorial/EquationsClassification.md",
	"tutorial/Membrane.md",
	"tutorial/HeatExchanger.md",
	"tutorial/Acoustics.md",
	"tutorial/ThermalConduction.md",
	"tutorial/FanBlade.md",
	"tutorial/RotatingHill.md",
	"tutorial/Elasticity.md",
	"tutorial/Stokes.md",
	"tutorial/NavierStokesProjection.md",
	"tutorial/NavierStokesNewton.md",
	"tutorial/ALargeFluidProblem.md",
	"tutorial/ComplexNumbers.md",
	"tutorial/OptimalControl.md",
	"tutorial/FlowWithShocks.md",
	"tutorial/HeatEquationOptimization.md",
	"tutorial/TimeDependentStokes.md",
	"tutorial/WifiPropagation.md"
	]
tutorialNamesList = [
	"Home",
	"Poisson",
	"Classification of the equations",
	"Membrane",
	"Heat exchanger",
	"Acoustics",
	"Thermal Conduction",
	"Irrotational Fan Blade Flow and Thermal effects",
	"Pure convection, The rotating hill",
	"The system of elasticity",
	"The system of Stokes for fluids",
	"A projection Algorithm for the Navier-Stokes equations",
	"Newton method for the steady Navier-Stokes equations",
	"A large fluid problem",
	"An example with complex numbers",
	"Optimal control",
	"A flow with shocks",
	"Time dependant schema optimization for heat equations",
	"A transient Stokes solver in matrix form",
	"Wifi Propagation"
	]
tutorialTODOFile = "tutorial/TODO.md"

## progressbar func
def writeProgressBar(file, Progress):
	if Progress == 0:
		file.write("<div class=\"progress progress-0\">\n")
	elif Progress < 20:
		file.write("<div class=\"progress progress-0plus\">\n")
	elif Progress < 40:
		file.write("<div class=\"progress progress-20plus\">\n")
	elif Progress < 60:
		file.write("<div class=\"progress progress-40plus\">\n")
	elif Progress < 80:
		file.write("<div class=\"progress progress-60plus\">\n")
	elif Progress < 100:
		file.write("<div class=\"progress progress-80plus\">\n")
	else:
		file.write("<div class=\"progress progress-100plus\">\n")

	file.write("\t<div class=\"progress-bar\" style=\"width:"+str(Progress)+"%\">\n")
	file.write("\t</div>\n")
	file.write("\t<span class=\"progress-label\">"+str(Progress)+"</span>\n")
	file.write("</div>\n\n")

## check TODO by lines
def checkTODOByLines(_FilesList, _NamesList, _TODOFile):
	TODOFile = open(docDirectory + _TODOFile, "w")
	TODOFile.write("<!--- THIS FILE IS AUTOMATICALY GENERATED --->\n")
	TODOFile.write("<!--- DO NOT EDIT --->\n\n")
	TODOFile.write("# TODO\n\n")
	for i in range(0, len(_FilesList)):
		file = open(docDirectory + _FilesList[i], "r")
		if file:
			TODOFile.write("## "+_NamesList[i]+"\n\n")

			EmptyLine = 0
			CurrentLine = 0
			NumberOfTODO = 0
			NumberOfError = 0
			TODOs = []
			for line in file:
				line = line.replace("\n", "")
				if not line:
					EmptyLine += 1
					CurrentLine += 1
				else:
					CurrentLine += 1
					if "$\\codered$" in line:
						TODOs.append("- [ ] line "+str(CurrentLine)+"\n")
						NumberOfTODO += 1
					elif "$\\codeerror$" in line:
						TODOs.append("- [ ] line "+str(CurrentLine)+"\n")
						TODOs.append("\n<span style=\"color:red; font-size:1.5em;\">This is a CodeError!</span>\n")
						print "A CodeError is present: file "+_FilesList[i]
						NumberOfError += 1

			TODOFile.write("Progression:\n")
			if NumberOfError != 0:
				Progression = 0
			else:
				Progression = 100
				if (CurrentLine-EmptyLine) > 0:
					Progression = 100 - (float(NumberOfTODO)/float(CurrentLine-EmptyLine))*100

			Progression = int(Progression)
			writeProgressBar(TODOFile, Progression)

			for TODO in TODOs:
				TODOFile.write(TODO)
			TODOFile.write("\n")

			file.close()

	TODOFile.close()

## check $\codered$
#index
checkTODOByLines(indexFilesList, indexNamesList, indexTODOFile)

#documentation
checkTODOByLines(documentationFilesList, documentationNamesList, documentationTODOFile)

#examples
checkTODOByLines(examplesFilesList, examplesNamesList, examplesTODOFile)

#introduction
checkTODOByLines(introductionFilesList, introductionNamesList, introductionTODOFile)
TODOFile = open(docDirectory + introductionTODOFile, "a")
TODOFile.write("\n## Documentation\n");
TODOFile.write("[TODO](../documentation/TODO)\n")
TODOFile.write("\n## Language references\n");
TODOFile.write("[TODO](../reference/TODO)\n")
TODOFile.write("\n## Tutorials\n");
TODOFile.write("[TODO](../tutorial/TODO)\n")
TODOFile.write("\n## Examples\n");
TODOFile.write("[TODO](../examples/TODO)\n")
TODOFile.close();

#reference (special TODO)
TODOFile = open(docDirectory + referenceTODOFile, "w")
TODOFile.write("<!--- THIS FILE IS AUTOMATICALY GENERATED --->")
TODOFile.write("<!--- DO NOT EDIT --->")
TODOFile.write("# TODO\n\n")
TODOFile.write("## Add [optional] tag to parameters.\n\n")
TODOFile.write("Progression:\n")
TODOFile.write("Unknown\n\n")
TODOFile.write("See [int2d](functions/#int2d) for example.\n\n")
for i in range(0, len(referenceFilesList)):
	file = open(docDirectory + referenceFilesList[i], "r")
	if file:
		TODOFile.write("## "+referenceNamesList[i]+"\n\n")

		CurrentItem = ""
		NumberOfItems = 0
		NumberOfTODO = 0
		TODOs = []
		for line in file:
			if "##" in line:
				CurrentItem = line
				NumberOfItems += 1
			if "$\\codered$" in line:
				Item = CurrentItem.replace("#", "")
				Item = Item.replace("\n", "")
				if Item[0] == " ":
					Item = Item[1:]
				Link = referenceFilesList[i].replace(".md", "")
				Link = Link.replace("reference/", "")
				TODOs.append("- [ ] ["+Item+"]("+Link+"/#"+Item+")\n")
				NumberOfTODO += 1

		TODOFile.write("Progression:\n")
		Progression = 100
		if NumberOfItems != 0:
			Progression = 100. - (float(NumberOfTODO)/float(NumberOfItems))*100.
		Progression = int(Progression)
		writeProgressBar(TODOFile, Progression)

		for TODO in TODOs:
			TODOFile.write(TODO)
		TODOFile.write("\n")

		file.close()

TODOFile.close()

#tutorial
checkTODOByLines(tutorialFilesList, tutorialNamesList, tutorialTODOFile)
