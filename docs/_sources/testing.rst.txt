.. _testing:

Testing
=======

Python testing
--------------

Testing of python scripts uses `pytest`_.

Set the :code:`PYTHONPATH` environment variable to the :code:`bin` directory where the scripts are stored, and then run the following command:

::

    python -m pytest -q tests/test_class.py

.. _pytest: https://docs.pytest.org/en/7.1.x/