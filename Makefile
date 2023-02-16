# options
.ONESHELL:


.PHONY: deps-install
deps-install:  ## install dependencies
	pip install poetry
	poetry install 
	poetry shell
