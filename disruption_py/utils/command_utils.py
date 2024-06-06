#!/usr/bin/env python3

import subprocess


def get_commit_hash():
	# setup commit hash
	try:
		commit_hash = subprocess.check_output(
            ['git', 'rev-parse', '--short', 'HEAD']
        ).decode('ascii').strip()
	except Exception as e:
		commit_hash = None
	return commit_hash
