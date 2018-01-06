======================
 Developing seqmagick
======================

Requirements
============

Note that building docs, publishing to pypi, etc require some
additional dependencies. It's best to work in a virtualenv::

  python3 -m venv py3-env
  source py3-env/bin/activate
  pip install -r requirements.txt

Git workflow
============

We aspire to more or less use a `feature branch workflow
<https://www.atlassian.com/git/workflows#!workflow-feature-branch>`_
for development. Briefly (for those working on the main fork):

* Features or bugfixes should start life as a GitHub issue
* Work on the feature occurs in a "feature branch" named like
  '%i-brief-description' % issue_number
* When completed (with tests passing) the feature branch is merged
  into dev (a pull request at this point might be appropriate if you
  want to request a code review).
* When it's time for a release, dev is merged into master (as a
  result, the head of the `master` branch is always on a release
  version).

versioning
==========

The package version is defined by the git tag (using
``git describe --tags --dirty``), in the form '<major>.<minor>.<bugfix>',
eg::

  git tag -a -m 'version 0.7.0' 0.7.0

Because setup.py determines the package version at the time the
package tarball is created, the repo must be clean (ie, no uncommitted
changes to versioned files) and there must be no further commits after
adding the tag when preparing to upload a tarball to PyPi.

preparing a release
===================

First, make sure you have committed all changes.

Run tests, and make sure docs build without errors::

  nosetests
  (cd docs && make html)

Push one last time to master to trigger tests on travis::

  git push origin master

Go to travis (https://travis-ci.org/fhcrc/seqmagick) and make sure the
tests have completed.

Add a new tag (see above). Push the tag to GitHub::

  git push --tags

Build and upload a tarball to PyPi::

  python setup.py clean
  rm -r build dist
  python setup.py sdist
  twine upload dist/*

Build and push docs to GutHub pages::

  (cd docs && make html)
  ghp-import --no-jekyll -p docs/_build/html
