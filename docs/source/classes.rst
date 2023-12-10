=======
Classes
=======

The classes are the classes are stored under the namespace ``curv3d``, and hence
is omitted from now on.

.. note::

   Most functions should also available in the Python module, and will follow the
   same name and parameters unless otherwise specified.

We can sort the classes by file:

#. :ref:`Utilities`
#. :ref:`Primitive Related`
#. :ref:`Object Related`
#. :ref:`Scene Related`


Utilities
=========
These define commonly used classes and types.


Base Types
----------

.. doxygentypedef:: real

.. doxygentypedef:: uint

Tuple Types
-----------

.. doxygentypedef:: curv3d::SPoint

.. doxygentypedef:: curv3d::SVector

.. doxygentypedef:: curv3d::Order

.. doxygentypedef:: curv3d::TwoIndex

.. doxygentypedef:: curv3d::GPoint

.. doxygentypedef:: curv3d::LPoint

.. doxygentypedef:: curv3d::Vector

.. doxygentypedef:: curv3d::Angles

Both these vector structs support all expected operator overloads.

TwoVector<T>
------------
.. doxygenstruct:: curv3d::TwoVector
    :members:
    :allow-dot-graphs:

ThreeVector<T>
--------------
.. doxygenstruct:: curv3d::ThreeVector
    :members:
    :allow-dot-graphs:


Primitive Related
=================

Surface
-------

.. doxygenclass:: curv3d::Surface
    :members:
    :private-members:
    :allow-dot-graphs:


Object Related
==============

.. doxygentypedef:: curv3d::Basis

.. doxygentypedef:: curv3d::PointMatrix

.. doxygentypedef:: curv3d::ScreenBasis

Object
------

.. doxygenclass:: curv3d::Object
    :members:
    :protected-members:
    :allow-dot-graphs:

Camera
------

.. doxygenclass:: curv3d::Camera
    :members:
    :protected-members:
    :allow-dot-graphs:

Scene Related
=============

Pixel
-----

.. doxygenstruct:: curv3d::Pixel
    :members:
    :allow-dot-graphs:

ImageBuffer
-----------

.. doxygenstruct:: curv3d::ImageBuffer
    :members:
    :allow-dot-graphs:

Scene
-----

.. doxygenclass:: curv3d::Scene
    :members:
    :private-members:
    :allow-dot-graphs:
