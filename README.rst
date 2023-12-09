======
curv3d
======
This repository is for `curv3d`, a Bézier surface based 3D engine.

Building
========
All that is required is a working Python installation. First, clone the repository

.. code-block:: sh

    git clone https://github.com/monokuchi/curv3d.git
    cd curv3d
    git submodule update --init --recursive

You will need `cmake` to continue. Next, we build the dependencies by running

.. code-block:: sh

    ./SETUP

Lastly, you can build the Python module by running the Makefile

.. code-block:: sh

    make python

Documentation
=============
To build the documentation, you simply only need to run the Makefile

.. code-block:: sh

    make docs
