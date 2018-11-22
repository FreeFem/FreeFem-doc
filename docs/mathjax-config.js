MathJax.Hub.Config({
	tex2jax: {
	inlineMath: [
		['$', '$'],
		["\\(", "\\)"]
	],
	processEscapes: true
	},
	jax: ["input/TeX", "output/SVG"],
	TeX: {
		TagSide: "right",
		TagIndent: ".8em",
		equationNumbers: {
			autoNumber: "AMS"
		},
		Macros: {
			codeerror: '{\\color{red}{ERROR}}',
			codered: '{\\color{red}{TODO}}',
			R: '{\\mathbb R}',
			N: '{\\mathbb N}',
			C: '{\\mathbb C}',
			P: '{\\mathbb P}',
			p: '{\\partial}',
			d: '{\\text d}',
			n: '{\\nabla}',
			T: '{\\mathbb T}',
			vecttwo: ['\\left|\\begin{array}{c}{#1}\\\\#2\\end{array}\\right.', 2],
			vectthree: ['\\left|\\begin{array}{c}{#1}\\\\#2\\\\#3\\end{array}\\right.', 3],
			MyInt: ['\\displaystyle{\\int_{#1} #2 }', 2]
		}
	}
});
