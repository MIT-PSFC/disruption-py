
PYLINT_DIRS := disruption_py examples tests

# git #

.PHONY: status fetch

status: fetch
	git status

fetch:
	git fetch -p -a

# poetry #

.PHONY: install uninstall lock update show

install:
	poetry install --with dev

uninstall:
	poetry env list | cut -d' ' -f1 | xargs poetry env remove

lock:
	poetry lock --no-update
	bash .pip/lock.sh
	git status

update:
	poetry update

show:
	poetry show --latest --why

# test #

.PHONY: quick test test-fast

quick:
	poetry run pytest -v tests/test_quick.py

test:
	poetry run pytest -v tests

test-fast:
	GITHUB_ACTIONS=1 poetry run pytest -v tests

# lint #

.PHONY: lint black pylint shellcheck yamllint

lint: black pylint shellcheck yamllint

black:
	poetry run black --version
	poetry run black --check .

pylint:
	poetry run pylint --version
	poetry run pylint $(PYLINT_DIRS)

pylint-only:
	poetry run pylint --disable=all --enable=$(CODE) $(PYLINT_DIRS)

shellcheck:
	poetry run shellcheck --version
	find -type f -not -path '*/.git/*' \
	| xargs grep -l '^#!/bin/bash' \
	| while read -r F; \
	do \
	   echo "--> $$F"; \
	   poetry run shellcheck "$$F"; \
	done

yamllint:
	poetry run yamllint --version
	find -type f -iname '*.yml' -or -iname '*.yaml' \
	| while read -r F; \
	do \
	   echo "--> $$F"; \
	   poetry run yamllint "$$F"; \
	done
