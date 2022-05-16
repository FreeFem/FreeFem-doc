# Sphinx installation for FreeFEM-doc

```
pip install -U Sphinx
```

## Extensions

```
pip install sphinxcontrib-inlinesyntaxhighlight sphinx-sitemap sphinxcontrib-svg2pdfconverter[CairoSVG]
```

## Usage

On root directory, use:
- `make htmlonly` to build the html website
- `make html` to build the html website and the [lunr](https://lunrjs.com/) index (search engine). [Nodejs](https://nodejs.org/en/) is required
- `make latex` to build the latex
