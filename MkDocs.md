# Install MkDocs with dependencies

## On Linux / MacOS

Install python:
```
apt-get install python
```

Install pip:
```
apt-get install python-pip
```

Update pip:
```
pip install --upgrade pip
```

Install MkDocs:
```
pip install mkdocs
```

Install Material theme for MkDocs:
```
pip install mkdocs-material
```

Install math extension:
```
pip install python-markdown-math
```

That's it! Use `mkdocs serve` for the development server or `mkdocs build` to build the static web pages.

## On Windows

Download latest release from [Python website](https://www.python.org/) and install it. ***Verify that environment variables export is enabled***.

Install MkDocs, Material and math extension:
```
python -m pip install mkdocs
python -m pip install mkdocs-material
python -m pip install python-markdown-math
```
