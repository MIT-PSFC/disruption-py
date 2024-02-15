import subprocess


def get_commit_hash():
	# setup commit hash
	try:
		commit_hash = subprocess.check_output(
			["git", "describe", "--always"],
			stdout=subprocess.DEVNULL,
			stderr=subprocess.STDOUT).strip()
	except Exception as e:
		# self.logger.warning("Git commit not found")
		commit_hash = 'Unknown'