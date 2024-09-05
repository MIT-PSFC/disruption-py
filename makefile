
PYLINT_DIRS := disruption_py examples tests
DELETE_OBJS := __pycache__ .pytest_cache

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

.PHONY: install uninstall reinstall lock update show

install:
	poetry install --with dev

uninstall:
	poetry env list | cut -d' ' -f1 | xargs poetry env remove

reinstall: uninstall install

lock:
	poetry lock --no-update
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

.PHONY: lint black isort pylint shellcheck yamllint

lint: black pylint shellcheck yamllint

black:
	poetry run black --version
	poetry run black --check .

isort:
	poetry run isort --version
	poetry run isort --check --profile black .

pylint:
	poetry run pylint --version
	find $(PYLINT_DIRS) -type f -name '*.py' \
	| xargs poetry run pylint

pylint-only:
	find $(PYLINT_DIRS) -type f -name '*.py' \
	| xargs poetry run pylint --disable=all --enable=$(CODE)

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
