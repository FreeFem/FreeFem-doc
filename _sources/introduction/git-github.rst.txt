Git & Github usage
==================

FreeFEM sources are publicly available on https://github.com/FreeFem/FreeFem-sources.

In order to contribute, you need to know how to use git (:code:`add`, :code:`commit`, :code:`push`) and Github (Fork, Pull Requests).

The FreeFEM source code is organized in branches:

    - :code:`master`. The master branch, represent the current stable version, used to build a new release

    - :code:`develop`. The developement branch, where all modifications take place
        Should be almost always usable

    - features branches, where specific long-term developments take place
        Do not use one of this branch

Contribution timeline
---------------------

    - Create a fork of the FreeFem-sources repository on your Github account
        Doc: |Fork_documentation|

        Direct fork link: |Fork_link|

    - Clone the fork (the FreeFem-sources repository on your account) on your computer.
        Change the branch to :code:`develop`
        
        :code:`git checkout develop`

    - Modify the code and use git commands to push your modifications to the fork, i.e.:
        :code:`git add somefile.cpp`

        :code:`git commit -m"my modification"`

        :code:`git push`

        Please, provide commit descriptions correctly describe your modifications

    - Create a pull request on FreeFem/FreeFem-sources, describing your modifications
        Doc: |PR_documentation|

.. warning:: All code modifications, even in a pull request, must be done in the `develop` branch

.. note:: Please make sure your code modification is well writen and formatted (you can use clang-format if necessary)

.. |Fork_documentation| raw:: html

    <a href="https://docs.github.com/en/get-started/quickstart/fork-a-repo" target="_blank">Github Fork documentation</a>

.. |Fork_link| raw:: html

    <a href="https://github.com/FreeFem/FreeFem-sources/fork" target="_blank">Fork link</a>

.. |PR_documentation| raw:: html

    <a href="https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request" target="_blank">Pull Request documentation</a>
