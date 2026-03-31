BRANCHES := alpha stable
.PHONY: $(BRANCHES) all

$(BRANCHES):
	cd $@/ && git fetch -p && git status && git log -1

$(BRANCHES:%=%-build):
	cd $(patsubst %-build,%,$@) && \
	git fetch  && \
	BRANCH=$$(git branch --show-current) && \
	LATEST=$$(git tag -l "$$BRANCH*" --format='%(refname:short)' --sort=-creatordate | head -n1) && \
	echo "$$BRANCH @ $$LATEST" && \
	git reset --hard "$$LATEST" && \
	git status && \
	git log -1 && \
	git clean -fdx && \
	./bootstrap && \
	cd python/MDSplus && \
	find /usr/local/mdsplus /fusion/usc/c8/opt/mdsplus/{alpha,stable} -path '*/MDSplus/_version.py' 2> /dev/null | \
	xargs ls -t | \
	head -n1 | \
	xargs ln -sfv && \
	realpath _version.py && \
	cat _version.py

all: $(BRANCHES:%=%-build)
