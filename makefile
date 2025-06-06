
# parameters #

PYLINT_DIRS := disruption_py examples tests
DELETE_OBJS := __pycache__ .pytest_cache

# environment #

GITHUB_ACTIONS ?= 0
ifeq ($(GITHUB_ACTIONS), true)
    CHECK_ARG := --check
    FORMAT_ARG := --output-format=github
endif

# git #

.PHONY: status fetch

status: fetch
	git status

fetch:
	git fetch -p -a

# clean #

.PHONY: clean-list clean-delete

clean-list:
	echo $(DELETE_OBJS) | xargs -n1 find -name
	find -type d -empty

clean-delete:
	echo $(DELETE_OBJS) | xargs -n1 find -name | xargs rm -rfv
	find -type d -empty -delete

# poetry #

.PHONY: install install-all uninstall reinstall lock update show

install:
	poetry install --with dev

install-all:
	poetry install --with dev,docs,lab

uninstall:
	poetry env list | cut -d' ' -f1 | xargs poetry env remove

reinstall: uninstall install

lock:
	poetry lock
	uv lock || true
	git status

update:
	poetry update
	uv lock --upgrade || true
	git status

show:
	poetry show --latest --why --top-level --with dev,lab,docs

# test #

.PHONY: quick test test-fast

quick:
	poetry run pytest -v tests/test_quick.py

test:
	poetry run pytest -v --durations=0 tests

test-fast:
	GITHUB_ACTIONS=true poetry run pytest -v tests

# lint #

.PHONY: lint check isort black ruff pylint pylint-only pylint-todos shellcheck yamllint toml-sort

lint: isort black ruff pylint shellcheck yamllint toml-sort

check:
	make lint GITHUB_ACTIONS=true

black:
	@[ "$(GITHUB_ACTIONS)" != "true" ] || \
	poetry run black --version
	poetry run black $(CHECK_ARG) .

isort:
	@[ "$(GITHUB_ACTIONS)" != "true" ] || \
	poetry run isort --version
	poetry run isort $(CHECK_ARG) --profile black .

pylint:
	@[ "$(GITHUB_ACTIONS)" != "true" ] || \
	poetry run pylint --version
	find $(PYLINT_DIRS) -type f -name '*.py' -not -empty \
	| xargs poetry run pylint -v --jobs 4 $(FORMAT_ARG)

pylint-only:
	find $(PYLINT_DIRS) -type f -name '*.py' -not -empty  \
	| xargs poetry run pylint -v --disable=all --enable=$(CODE)

pylint-todos:
	CODE=fixme make pylint-only

ruff:
	@[ "$(GITHUB_ACTIONS)" != "true" ] || \
	poetry run ruff --version
	find $(PYLINT_DIRS) -type f -name '*.py' -not -empty \
	| xargs poetry run ruff check $(FORMAT_ARG)

shellcheck:
	@[ "$(GITHUB_ACTIONS)" != "true" ] || \
	shellcheck --version
	find -type f -not -path '*/.git/*' -not -path '*/.venv/*' \
	| xargs grep -l '^#!/bin/bash' \
	| while read -r F; \
	do \
	   echo "--> $$F"; \
	   shellcheck "$$F"; \
	done

yamllint:
	@[ "$(GITHUB_ACTIONS)" != "true" ] || \
	poetry run yamllint --version
	find -type f -iname '*.y*ml' -not -empty -not -path '*/.venv/*' \
	| while read -r F; \
	do \
	   echo "--> $$F"; \
	   poetry run yamllint "$$F"; \
	done

toml-sort:
	@[ "$(GITHUB_ACTIONS)" != "true" ] || \
	poetry run toml-sort --version
	find -maxdepth 1 -type f -iname '*.toml' -not -empty -not -path '*/.venv/*'\
	| while read -r F; \
	do \
	   echo "--> $$F"; \
	   poetry run toml-sort $(CHECK_ARG) "$$F"; \
	done
	find -mindepth 2 -type f -iname '*.toml' -not -empty -not -path '*/.venv/*' \
	| while read -r F; \
	do \
	   echo "--> $$F"; \
	   poetry run toml-sort $(CHECK_ARG) --all "$$F"; \
	done