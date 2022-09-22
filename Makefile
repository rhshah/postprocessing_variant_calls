# options
.ONESHELL:


.PHONY: deps-install
deps-install:  ## install dependencies
	pip install poetry
	poetry config virtualenvs.create false
	poetry install 

.PHONY: deps-update
deps-update:
	pip install poetry
	poetry config virtualenvs.create false
	poetry lock 
	poetry export --format requirements.txt --output requirements.txt --without-hashes

requirements.txt: poetry.lock
	poetry export --format requirements.txt --output requirements.txt --without-hashes

requirements-dev.txt: poetry.lock
	poetry export --dev --format requirements.txt --output requirements-dev.txt --without-hashes
