.. role:: freefem(code)
  :language: freefem

I/O
===

See :ref:`I/O example <exampleIO>`

See :ref:`File stream example <exampleFileStream>`.

cout
----

Standard C++ output device (default: console).

.. code-block:: freefem
   :linenos:

   cout << "Some text" << endl;

cin
---

Standard C++ input device (default: keyboard).

.. code-block:: freefem
   :linenos:

   cin >> var;

endl
----

End of line.

.. code-block:: freefem
   :linenos:

   cout << "Some text" << endl;

ifstream
--------

Open a file in read mode.

.. code-block:: freefem
   :linenos:

   ifstream file("file.txt");

.. note:: A file is closed at the end of a block.

ofstream
--------

Open a file in write mode.

.. code-block:: freefem
   :linenos:

   ofstream file("file.txt");

.. note:: A file is closed at the end of a block.

append
------

Append data to an existing file.

.. code-block:: freefem
   :linenos:

   ofstream file("file.txt", append);

binary
------

Write a file in binary.

.. code-block:: freefem
   :linenos:

   ofstream file("file.btxt", binary);

seekg
-----

Set the file position.

.. code-block:: freefem
   :linenos:

   file.seekg(Pos);

tellg
-----

Get the file position.

.. code-block:: freefem
   :linenos:

   int Pos = file.tellg();

flush
-----

Flush the buffer of the file.

.. code-block:: freefem
   :linenos:

   file.flush

getline
-------

Get the current line.

.. code-block:: freefem
   :linenos:

   string s;
   getline(file, s);

Output format
-------------

In the descriptions below, ``f`` is an output stream, for example :freefem:`cout` or a :freefem:`ofstream`.

All this methods, excepted the first, return a stream, so they can be chained:

.. code-block:: freefem
   :linenos:

   cout.scientific.showpos << 3 << endl;

precision
~~~~~~~~~

Set the number of digits printed to the right of the decimal point.
This applies to all subsequent floating point numbers written to that output stream.
However, this won’t make floating-point “integers" print with a decimal point.
It’s necessary to use :freefem:`fixed` for that effect.

.. code-block:: freefem
   :linenos:

   int np = f.precision(n)

scientific
~~~~~~~~~~

Formats floating-point numbers in scientific notation

.. code-block:: freefem
   :linenos:

   f.scientific

fixed
~~~~~

Used fixed point notation for floating-point numbers.
Opposite of scientific.

.. code-block:: freefem
   :linenos:

   f.fixed

showbase
~~~~~~~~

Converts insertions to an external form that can be read according to the ``C++`` lexical conventions for integral constants.
By default, showbase is not set.

.. code-block:: freefem
   :linenos:

   f.showbase

noshowbase
~~~~~~~~~~

Unset :freefem:`showbase` flags.

.. code-block:: freefem
   :linenos:

   f.noshowbase

showpos
~~~~~~~

Inserts a plus sign (+) into a decimal conversion of a positive integral value.

.. code-block:: freefem
   :linenos:

   f.showpos

noshowpos
~~~~~~~~~

Unset :freefem:`showpos` flags.

.. code-block:: freefem
   :linenos:

   f.noshowpos

default
~~~~~~~

Reset all the previous flags to the default expect precision.

.. code-block:: freefem
   :linenos:

   f.default

setw
~~~~

Behaves as if member width were called with ``n`` as argument on the stream on which it is inserted as a manipulator (it can be inserted on output streams).

.. code-block:: freefem
   :linenos:

   f.setw(n)
