#!/bin/sh

VERSION=2.0
make html && scp -r _build/html admol.org:/home/www/xcfun/$VERSION
