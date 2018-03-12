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
        	R: '{\\mathbb R}',
			P: '{\\mathbb P}',
        	p: '{\\partial}',
			d: '{\\text d}',
			n: '{\\nabla}',
			T: '{\\mathbb T}'
    	}
    }
});
