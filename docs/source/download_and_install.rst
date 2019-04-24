
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

Use the docker images
#########################
The pre-built docker image is easy to try out.
The docker images can be found `here <https://cloud.docker.com/repository/docker/zhouwanding/transvar/general>`__

Assuming the existence of `~/references/hg38/hg38.fa` and
`~/references/hg38/hg38.fa.fai`. 

Without downloading anything, the transvar docker has pre-built hg38
annotation.
.. code:: bash

	docker run -v ~/references/hg38:/ref -ti zhouwanding/transvar:latest transvar panno -i PIK3CA:p.E545K --ensembl --reference /ref/hg38.fa

To use other genome build, one needs to download annotations. Here I
am using `~/test` as an example of local path for storing the transvar
annotations. Note that this local path needs be imaged to `/anno`
inside the docker image. This is done by (showing hg19)
.. code:: bash
          
  docker run -v ~/test:/anno -ti zhouwanding/transvar:latest transvar config --download_anno --refversion hg19 --skip_reference

Now one can use hg19, but note again one needs to image the path of
downloaded annotation to `/anno`. One also needs the fa-indexed
reference.
.. code:: bash
          
  docker run -v ~/test:/anno -v ~/references/hg19:/ref -ti zhouwanding/transvar:latest transvar panno -i PIK3CA:p.E545K --ensembl --reference /ref/hg19.fa


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

