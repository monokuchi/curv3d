=======
Classes
=======

There are a few classes used for storing mesh data. All of the classes are stored under
the namespace ``Geometry``, and hence is ommitted from now on.

.. note::

   Most functions are also available in the Python module, however, the ``camelCase``
   naming is replaced with ``snake_case``. For example, in C++, one would write
   ``mesh.exportFacesSTL(...)``, but in Python, ``mesh.export_faces_stl(...)`` would be
   correct.

There are two main components:

#. :ref:`Mesh Related`
#. :ref:`IO Related`

Mesh Related
============
These classes relate primarily to the actual geometry mesh data itself.

ThreeVector<T>
--------------
.. doxygenclass:: Geometry::ThreeVector
    :members:
    :allow-dot-graphs:

.. doxygentypedef:: Geometry::Vector
.. doxygentypedef:: Geometry::Indices

Face
----
.. doxygenclass:: Geometry::Face
    :members:

Element
-------
.. doxygenclass:: Geometry::Element
    :members:


BoundaryCondition
-----------------
.. doxygenenum:: Geometry::BoundaryCondition

VolumeMeshInfo
--------------
.. doxygenstruct:: Geometry::VolumeMeshInfo
    :members:

VolumeMesh
----------
.. doxygenclass:: Geometry::VolumeMesh
    :members:
    :allow-dot-graphs:


IO Related
==========

EXODUS
------
.. doxygenclass:: Geometry::IO::EXODUS
    :members:

JSON
----

.. doxygenstruct:: Geometry::IO::JSON::ElementSet
    :members:

.. doxygenstruct:: Geometry::IO::JSON::FaceSet
    :members:

.. doxygenstruct:: Geometry::IO::JSON::HybridMesh
    :members:

.. doxygenstruct:: Geometry::IO::JSON::Config
    :members:
