.. role:: freefem(code)
  :language: freefem

Operators
=========

Addition operator +
-------------------

.. code-block:: freefem
   :linenos:

   real a = 1. + 2.;

Works for :freefem:`int`, :freefem:`real`, :freefem:`complex`, :freefem:`string`, :freefem:`mesh`, :freefem:`mesh3`, array.

Increment operator ++
---------------------

Pre-increment:

.. code-block:: freefem
   :linenos:

   int i = 0;
   ++i;

Post-increment:
.. code-block:: freefem
   :linenos:

   int i = 0;
   i++;

Substraction operator -
-----------------------

.. code-block:: freefem
   :linenos:

   real a = 1. - 2.;

Works for :freefem:`int`, :freefem:`real`, :freefem:`complex`, array.

Decrement operator --
---------------------

Pre-decrement:
.. code-block:: freefem
   :linenos:

   int i = 0;
   --i;

Post-decrement:
.. code-block:: freefem
   :linenos:

   int i = 0;
   i--;

Multiplication operator *
-------------------------

.. code-block:: freefem
   :linenos:

   real[int] b;
   matrix A
   real[int] x = A^-1*b;

Works for :freefem:`int`, :freefem:`real`, :freefem:`complex`, array, :freefem:`matrix`.

Equal operator =
----------------

.. code-block:: freefem
   :linenos:

   real a = 1.;

Comparison operator ==
----------------------

.. code-block:: freefem
   :linenos:

   real a = 1.;
   real b = 1.;

   cout << (a == b) << endl;

Comparison operator !=
----------------------

.. code-block:: freefem
   :linenos:

   real a = 1.;
   real b = 2.;

   cout << (a != b) << endl;

Comparison operator <, <=
-------------------------

.. code-block:: freefem
   :linenos:

   real a = 1.;
   real b = 2.;

   cout << (a < b) << endl;
   cout << (a <= b) << endl;

Comparison operator >, >=
-------------------------

.. code-block:: freefem
   :linenos:

   real a = 3.;
   real b = 2.;

   cout << (a > b) << endl;
   cout << (a >= b) << endl;

Term by term multiplication .*
------------------------------

.. code-block:: freefem
   :linenos:

   matrix A = B .* C;

Division operator /
-------------------

.. code-block:: freefem
   :linenos:

   real a = 1. / 2.;

Works for :freefem:`int`, :freefem:`real`, :freefem:`complex`.

Term by term division ./
------------------------

.. code-block:: freefem
   :linenos:

   matrix A = B ./ C;

Remainder from the division %
-----------------------------

.. code-block:: freefem
   :linenos:

   int a = 1 % 2;

Works for :freefem:`int`, :freefem:`real`.

Power operator ^
----------------

.. code-block:: freefem
   :linenos:

   real a = 2.^2;

Works for :freefem:`int`, :freefem:`real`, :freefem:`complex`, :freefem:`matrix`.

Inverse of a matrix ^-1
-----------------------

.. code-block:: freefem
   :linenos:

   real[int] Res = A^-1 * b;

.. warning:: This operator can not be used to directly create a matrix, see :ref:`Matrix inversion <exampleMatrixInversion>`.

Transpose operator '
--------------------

.. code-block:: freefem
   :linenos:

   real[int] a = b';

Works for array and :freefem:`matrix`.

.. note:: For :freefem:`matrix<complex>`, the ::freefem`'` operator return the Hermitian tranpose.

Tensor scalar product :
-----------------------

.. math::


   A:B = \sum_{i,j}{A_{ij}B_{ij}}

C++ arithmetical if expression ? :
--------------------------------------

``a ? b : c`` is equal to ``b`` if the ``a`` is true, ``c`` otherwise.

.. tip:: Example with :freefem:`int`

   .. code-block:: freefem
      :linenos:

      int a = 12; int b = 5;

      cout << a << " + " << b << " = " << a + b << endl;
      cout << a << " - " << b << " = " << a - b << endl;
      cout << a << " * " << b << " = " << a * b << endl;
      cout << a << " / " << b << " = " << a / b << endl;
      cout << a << " % " << b << " = " << a % b << endl;
      cout << a << " ^ " << b << " = " << a ^ b << endl;
      cout << "( " << a << " < " << b << " ? " << a << " : " << b << ") = " << (a < b ? a : b) << endl;

   The output of this script is:

   .. code-block:: bash

      12 + 5 = 17
      12 - 5 = 7
      12 * 5 = 60
      12 / 5 = 2
      12 % 5 = 2
      12 ^ 5 = 248832
      ( 12 < 5 ? 12 : 5) = 5

.. tip:: Example with :freefem:`real`

   .. code-block:: freefem
      :linenos:

      real a = qsrt(2.); real b = pi;

      cout << a << " + " << b << " = " << a + b << endl;
      cout << a << " - " << b << " = " << a - b << endl;
      cout << a << " * " << b << " = " << a * b << endl;
      cout << a << " / " << b << " = " << a / b << endl;
      cout << a << " % " << b << " = " << a % b << endl;
      cout << a << " ^ " << b << " = " << a ^ b << endl;
      cout << "( " << a << " < " << b << " ? " << a << " : " << b << ") = " << (a < b ? a : b) << endl;

   The output of this script is:

   .. code-block:: bash

      1.41421 + 3.14159 = 4.55581
      1.41421 - 3.14159 = -1.72738
      1.41421 * 3.14159 = 4.44288
      1.41421 / 3.14159 = 0.450158
      1.41421 % 3.14159 = 1
      1.41421 ^ 3.14159 = 2.97069
