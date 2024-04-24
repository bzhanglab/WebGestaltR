#!/bin/sh -e

# This script is used to update the vendored dependencies in the `vendor` directory.
# Code from github.com/yutannihilation/string2path/ under the MIT license.

cargo vendor

# c.f. https://reproducible-builds.org/docs/archives/
gtar \
  --sort=name \
  --mtime='1970-01-01 00:00:00Z' \
  --owner=0 \
  --group=0 \
  --numeric-owner \
  --xz \
  --create \
  --file=vendor.tar.xz \
  vendor

echo
echo
echo "#############################################"
echo "#                                           #"
echo "#  UPDATE src/config.toml !!!  #"
echo "#                                           #"
echo "#############################################"