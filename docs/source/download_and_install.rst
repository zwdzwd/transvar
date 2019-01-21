
********************
Download and Install
********************

Install using pip
###################

.. code:: bash
   
   sudo pip install transvar

   
or locally

.. code:: bash
          
   pip install --user transvar


to upgrade from a previous version

.. code:: bash

   pip install -U transvar

Try out the docker image
#########################
Also try out using the pre-built docker image
Assuming the existence of `~/references/hg38/hg38.fa` and `~/references/hg38/hg38.fa.fai`

.. code:: bash

	docker pull zhouwanding/transvar:2.4.6
	docker run -v ~/references/hg38:/data -ti zhouwanding/transvar:2.4.6 transvar panno -i PIK3CA:p.E545K --ensembl --reference /data/hg38.fa


Download the program
#######################

Current release
^^^^^^^^^^^^^^^^^

Latest release is available `here <https://github.com/zwdzwd/transvar/releases/latest>`__

For all previous versions, see `here <https://github.com/zwdzwd/transvar/releases>`__

Other old stable releases
^^^^^^^^^^^^^^^^^^^^^^^^^^

+ stable 2.0.x version `v2.0.12.20150626 <https://github.com/zwdzwd/transvar/archive/v2.0.12.20150626.zip>`__
+ stable 1.x version `v1.40 <https://github.com/zwdzwd/transvar/archive/v1.40.zip>`__

Dependency
############

The only requirement for building TransVar are Python 2.7 and a reasonably modern C compiler such as gcc.

Install from source
######################

Local install
^^^^^^^^^^^^^^^^

.. code:: bash

   python setup.py install --prefix [folder]

The installation will create two subfolders: ``[folder]/lib`` (which would contain libraries) and ``[folder]/bin`` (which would contain transvar executable).

When you run transvar, make sure ``[folder]/lib/python2.7/site-packages`` is in your PYTHONPATH. In some occasions, you need to ``mkdir -p [folder]/lib/python2.7/site-packages`` to make sure it exists before you could run ``setup.py``.
You can add it by putting

.. code:: bash

   export PYTHONPATH=$PYTHONPATH:[folder]/lib/python-2.7/site-packages/

to your ``.bashrc`` or ``.profile`` depending on your OS.

The installed executable is **[folder]/bin/transvar**.

System-wise install (need root)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: bash

   sudo python setup.py install

