#!/usr/bin/python
"""
Add or update the license information in the header of source files.

The script reads in the .gitattributes file, located in the project root
directory, to figure out which files need to be inspected and which license
header they need to have.
For example:

src/pedra/pedra_dlapack.F90 !licensefile
src/solver/*.hpp licensefile=.githooks/LICENSE-C++

The first line specifies that the file in src/pedra/pedra_dlapack.F90 should
not be touched, while the second line states that all .hpp files in src/solver
should get an header from the template in .githooks/LICENSE-C++
Location of files in .gitattributes are always specified with respect
to the project root directory.

The script reads in the appropriate license header template and prepares an
header with the correct year and authors information. These are encoded in the
variables YEAR and AUTHORS. The latter has to be modified by hand.
"""

from datetime import date
import glob
import os
import re
import shutil
import tempfile


def add_header(filepath, header, YEAR, AUTHORS):
    """
    Add or update header in source file
    """
    tmpdir = tempfile.gettempdir()
    tmpfil = os.path.join(tmpdir, os.path.basename(filepath) + '.bak')
    shutil.copy2(filepath, tmpfil)
    with open(tmpfil, 'r') as tmp:
        inpt = tmp.readlines()
        output = []

        # Check if header is already present
        present = re.compile(
            'XCFun, an arbitrary order exchange-correlation library')
        if list(filter(present.search, inpt)):
            # Check if year and authors in current file are up to date
            toupdate = re.compile(r'{0} (?!{1} {2}).*\n'.format(
                'Copyright \(C\)', YEAR, AUTHORS))
            if list(filter(toupdate.search, inpt)):
                print(('Updating header in {}'.format(filepath)))
                # Check to preserve '#!' at the top of the file
                if len(inpt) > 0 and inpt[0].startswith('#!'):
                    output.append(inpt[0] + '\n')
                    inpt = inpt[1:]
                regex = re.compile(r'Copyright \(C\).*\n')
                repl = r'Copyright (C) ' + YEAR + ' ' + AUTHORS + '\n'
                output.extend([re.sub(regex, repl, x) for x in inpt])
        else:
            print(('Adding header in {}'.format(filepath)))
            # Check to preserve '#!' at the top of the file
            if len(inpt) > 0 and inpt[0].startswith('#!'):
                output.append(inpt[0] + '\n')
                inpt = inpt[1:]
            output.append(header)
            for line in inpt:
                output.append(line)

        if output:
            try:
                f = open(filepath, 'w')
                f.writelines(output)
            except IOError as err:
                print(('Something went wrong trying to add header to {}: {}'.
                       format(filepath, err)))
            finally:
                f.close()
        os.remove(tmpfil)


def prepare_header(stub, YEAR, AUTHORS):
    """
    Update year and author information in license header template
    """
    with open(stub, 'r') as l:
        header = l.read()
        # Insert correct YEAR and AUTHORS in stub
        rep = {'YEAR': YEAR, 'AUTHORS': AUTHORS}
        rep = dict((re.escape(k), v) for k, v in rep.items())
        pattern = re.compile("|".join(list(rep.keys())))
        header = pattern.sub(lambda m: rep[re.escape(m.group(0))], header)
    return header


def file_license(attributes):
    """
    Obtain dictionary { file : license } from .gitattributes
    """
    file_license = {}
    with open(attributes, 'r') as f:
        # Read in .gitattributes
        tmp = f.read()
        # Removing all comment lines and other attributes
        pattern = re.compile(r'(?m)^\#.*\n?|^((?!licensefile).)*$')
        gitattributes = re.sub(pattern, '', tmp).split()
        # Obtain list of files
        fil = [x for x in gitattributes if not 'licensefile' in x]
        # Remove licensefile= from strings
        lic = [
            re.sub(r'licensefile\=', '', x) for x in gitattributes
            if 'licensefile' in x
        ]
        # Create list of blacklisted files
        blacklist = [
            fname for key, value in list(dict(list(zip(fil, lic))).items())
            if value == '!licensefile' for fname in glob.glob(key)
        ]
        # Now create a dictionary with the files to be considered for
        # license header manipulation
        file_license = {
            key: value
            for k, value in list(dict(list(zip(fil, lic))).items())
            for key in glob.glob(k) if key not in blacklist
        }
    return file_license


def license_maintainer():
    """
    Maintain license header in source files
    """
    YEAR = str(date.today().year)
    AUTHORS = 'Ulf Ekström and contributors.'

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root_dir = os.path.abspath(os.path.join(script_dir, os.pardir))

    headerize = file_license(os.path.join(project_root_dir, '.gitattributes'))

    for fname, license in list(headerize.items()):
        # Prepare header
        header = prepare_header(
            os.path.join(project_root_dir, license), YEAR, AUTHORS)
        add_header(fname, header, YEAR, AUTHORS)


if __name__ == '__main__':
    license_maintainer()
