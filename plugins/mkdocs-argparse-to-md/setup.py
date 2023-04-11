from setuptools import setup, find_packages

setup(
    name="mkdocs-argparse-to-md",
    version="0.1.0",
    description="A custom MkDocs plugin for auto-generating Markdown tables from argparse arguments",
    author="Herbert Turner",
    author_email="hmturner@mit.edu",
    packages=['argparse_to_md'],
    install_requires=["mkdocs"],
    entry_points={
        "mkdocs.plugins": [
            "argparse_to_md = argparse_to_md.plugin:ArgparseToMDPlugin"
        ]
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)
