.. role:: freefem(code)
  :language: freefem

Loops
=====

See :ref:`Loop example <exampleLoop>`.

.. _loopFor:

for
---

For loop.

.. code-block:: freefem
   :linenos:

   for (int i = 0; i < N; ++i){
       ...
   }

.. _loopIf:

if
--

If condition.

.. code-block:: freefem
   :linenos:

   if (condition){
       ...
   }
   else{
       ...
   }

else
----

See :ref:`if <loopIf>`.

while
-----

While loop.

.. code-block:: freefem
   :linenos:

   while (condition){
       ...
   }

continue
--------

Continue a loop.

.. code-block:: freefem
   :linenos:

   for (int i = 0; i < N; ++i){
       ...
       if (condition) continue;
       ...
   }

break
-----

Break a loop.

.. code-block:: freefem
   :linenos:

   while (condition1){
       ...
       if (condition) break;
       ...
   }

.. _loopTry:

try
---

Try a part of code.

.. code-block:: freefem
   :linenos:

   try{
       ...
   }
   catch(...){
       ...
   }

See :ref:`Basic error handling example <exampleBasicErrorHandling>` and :ref:`Error handling example <exampleErrorHandling>`.

catch
-----

Catch an error, see :ref:`try <loopTry>`

Implicit loop
-------------

Array with one index:

.. code-block:: freefem
   :linenos:

   for [i, ai : a]

If :freefem:`real[int] a(10)`, then ``i=0:9`` and ``ai`` is a reference to ``a[i]``.

Array with two indices or matrix:

.. code-block:: freefem
   :linenos:

   for [i, j, aij : a]

If :freefem:`real[int] a(10, 11)`, then ``i=0:9``, ``j=1:10`` and ``aij`` is a reference to ``a(i, j)``.

See :ref:`Implicit loop example <exampleImplicitLoop>`.
