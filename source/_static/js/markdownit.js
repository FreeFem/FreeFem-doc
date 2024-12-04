var defaults = {
	html: false, // Enable HTML tags in source
	xhtmlOut: false, // Use '/' to close single tags (<br />)
	breaks: false, // Convert '\n' in paragraphs into <br>
	langPrefix: 'language-', // CSS language prefix for fenced blocks
	linkify: true, // autoconvert URL-like texts to links
	typographer: true, // Enable smartypants and other sweet transforms
	_highlight: true,
	_strict: false,
	_view: 'html' // html / src / debug
};

defaults.highlight = function (str, lang) {
	var esc = md.utils.escapeHtml;
	try {
		var html = "<div class = \"highlight\"><pre class=\"cm-s-default\">";
		const outf = function(token, tokenClass) {
			html = html + "<span class=\"cm-" + tokenClass + "\">" + token + "</span>"
		}
		CodeMirror.runMode(str, 'text/x-ff++src',outf);
		html = html + "</pre></div>";
		return html;
	} catch (__) {}
};

var md = window.markdownit(defaults);
md.use(window.markdownitMathjax(), {options :{
	beforeMath: '\n',
	afterMath: '\n',
	beforeInlineMath: '$',
	afterInlineMath: '$',
	beforeDisplayMath: '$$',
	afterDisplayMath: '$$'}
});
