# options
.ONESHELL:


.PHONY: deps-install
deps-install:  ## install dependencies
	pip install poetry
	poetry install 

.PHONY: deps-update
deps-update:
	pip install poetry
	poetry lock 
	poetry export --format requirements.txt --output requirements.txt --without-hashes

requirements.txt: poetry.lock
	poetry export --format requirements.txt --output requirements.txt --without-hashes

requirements-dev.txt: poetry.lock
	poetry export --dev --format requirements.txt --output requirements-dev.txt --without-hashes

init-submodule:
	git submodule update --init

clean-submodule:
	git submodule deinit -f .

init-submodule:
	git submodule update --init

update-submodule:
	git submodule update --recursive --remote

init-tools:
	cd python_bed_lookup \
	python setup.py build \
	python setup.py install \
