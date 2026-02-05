#!/usr/bin/env -S uv run --script
# /// script
# dependencies = [
#   "netcdf4",
#   "requests",
#   "xarray",
# ]
# ///

"""
simple script to check that imas links in an xarray dataset are valid.
"""

import sys

import requests
import xarray as xr
from colorama import Fore

cache = {}

for k, da in xr.open_dataset(sys.argv[1]).items():

    print(k)

    imas = da.attrs.get("imas")
    url = da.attrs.get("url")
    if not imas or not url:
        continue

    base, anchor = url.split("#")
    found = f'href="#{anchor}"' in cache.get(base, requests.get(base).text)

    print(f"\timas = \x1b[{'32m' if found else '31m'}{imas}\x1b[39m")
    print(f"\turl  = {url}")
