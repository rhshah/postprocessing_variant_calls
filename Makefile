## options
# based on https://tech.davis-hansson.com/p/make/
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules
SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.ONESHELL:
.DEFAULT_GOAL := help
.DELETE_ON_ERROR:

## variables
POETRY_VERSION = 1.1.13
ENVIRONMENT ?= dev
ARGS =
APP_NAME = $(shell python -m src.config app_name)
SOURCE_DIR := src
TEST_DIR := tests

IMAGE_HOST = $(shell python -m src.config image_host)
IMAGE_REPO = $(shell python -m src.config image_repo)
IMAGE_NAME = $(IMAGE_HOST)/$(IMAGE_REPO)/$(APP_NAME)
IMAGE_TAG ?= latest

## formula

# based on https://marmelab.com/blog/2016/02/29/auto-documented-makefile.html
.PHONY: help
help:  ## print help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

## dependencies

.PHONY: deps-install
deps-install:  ## install dependencies
	pip install poetry==$(POETRY_VERSION)
	poetry config virtualenvs.create false
	poetry install --no-interaction 

.PHONY: deps-update
deps-update:
	pip install poetry==$(POETRY_VERSION)
	poetry config virtualenvs.create false
	poetry lock 
	poetry export --format requirements.txt --output requirements.txt --without-hashes

requirements.txt: poetry.lock
	poetry export --format requirements.txt --output requirements.txt --without-hashes

requirements-dev.txt: poetry.lock
	poetry export --dev --format requirements.txt --output requirements-dev.txt --without-hashes
