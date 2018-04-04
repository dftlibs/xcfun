#!/usr/bin/env bash

#
#  PCMSolver, an API for the Polarizable Continuum Model
#  Copyright (C) 2016 Roberto Di Remigio, Luca Frediani and collaborators
#
#  This file is part of PCMSolver.
#
#  PCMSolver is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  PCMSolver is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with PCMSolver.  If not, see <http://www.gnu.org/licenses/>.
#
#  For information on the complete list of contributors to the
#  PCMSolver API, see: <http://pcmsolver.readthedocs.io/>
#
# Provide the canonicalize filename (physical filename with out any symlinks)
# like the GNU version readlink with the -f option regardless of the version of
# readlink (GNU or BSD).

# This file is part of a set of unofficial pre-commit hooks available
# at github.
# Link:    https://github.com/githubbrowser/Pre-commit-hooks
# Contact: David Martin, david.martin.mailbox@googlemail.com

###########################################################
# There should be no need to change anything below this line.

# Canonicalize by recursively following every symlink in every component of the
# specified filename.  This should reproduce the results of the GNU version of
# readlink with the -f option.
#
# Reference: http://stackoverflow.com/questions/1055671/how-can-i-get-the-behavior-of-gnus-readlink-f-on-a-mac
canonicalize_filename () {
    local target_file="$1"
    local physical_directory=""
    local result=""

    # Need to restore the working directory after work.
    local working_dir="`pwd`"

    cd -- "$(dirname -- "$target_file")"
    target_file="$(basename -- "$target_file")"

    # Iterate down a (possible) chain of symlinks
    while [ -L "$target_file" ]
    do
        target_file="$(readlink -- "$target_file")"
        cd -- "$(dirname -- "$target_file")"
        target_file="$(basename -- "$target_file")"
    done

    # Compute the canonicalized name by finding the physical path
    # for the directory we're in and appending the target file.
    physical_directory="`pwd -P`"
    result="$physical_directory/$target_file"

    # restore the working directory after work.
    cd -- "$working_dir"

    echo "$result"
}
