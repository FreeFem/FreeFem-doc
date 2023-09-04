# Sphinx installation for FreeFEM-doc

```
pip install -U Sphinx
```

## Extensions

```
pip install sphinxcontrib-inlinesyntaxhighlight sphinx-sitemap
python3 -c "import site; print(site.getsitepackages()[0]+'/sphinxcontrib/inlinesyntaxhighlight.py')" | xargs sed -I '' -e 's/, BaseTranslator//g' -e "s/self.highlightlang/'python'/g"

```

## Usage

On root directory, use:
- `make htmlonly` to build the html website
- `make html` to build the html website and the [lunr](https://lunrjs.com/) index (search engine). [Nodejs](https://nodejs.org/en/) is required
- `make latex` to build the latex

## Install Node.js for the lunr index

In order to build the html website with the [lunr](https://lunrjs.com/) index, you need the [ Node.js](https://nodejs.org/) framework and `npm`.  

- Node installation on Linux

  ```
  sudo apt install nodejs
  sudo apt install npm
  ```
- Node installation on MacOS
  ```
  brew install node
  ```