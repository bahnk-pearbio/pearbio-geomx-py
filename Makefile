
PACKAGE=pearbio_geomx_py

.PHONY: .FORCE

all: .FORCE zip

install:
	python3 -m pip install --no-cache-dir .

zip: readme test clean
	cd ../ && zip -r $(PACKAGE).zip $(PACKAGE)

readme:
	pandoc --toc -o README.pdf README.md

clean:
	find $(PACKAGE) -name "__pycache__" -type d -exec rm -rfv {} \; || true
	rm -rfv build/
	rm -rfv dist/
	rm -rfv $(PACKAGE).egg-info/
	rm -rfv docs/build/
	rm -f docs/make.bat
	rm -rfv .pytest_cache

build:
	python3 -m build

doc:
	cd docs && make html

test: pylint
	pytest --verbose

pylint: black
	pylint setup.py
	pylint --recursive y $(PACKAGE)

black:
	black setup.py
	black $(PACKAGE)
