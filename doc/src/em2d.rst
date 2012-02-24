.. -*- coding: utf-8 -*-
.. $Id$

===================
2次元PICコード
===================
2次元コードはMPI並列版のみです。

.. blockdiag::

   diagram {
    node_width = 150;
    node_height = 25;
    span_width = 15;
    span_height = 15;

    "$PCANS_DIR/em2d_mpi/" -- "Makefile";
    "$PCANS_DIR/em2d_mpi/" -- "Makefile_inc";
    "$PCANS_DIR/em2d_mpi/" -- "common/"; 
    "$PCANS_DIR/em2d_mpi/" -- "moment/";
    "$PCANS_DIR/em2d_mpi/" -- "psd/";
    "$PCANS_DIR/em2d_mpi/" -- "md_wave/"
    "$PCANS_DIR/em2d_mpi/" -- "md_shock/" ;
    "$PCANS_DIR/em2d_mpi/" -- "md_kh/" ;
   }

基本的な使い方は :ref:`1次元シリアル版 <em1d>` と同じです。


課題例
======
.. toctree::
   :maxdepth: 2

   em2d_kh
.. em2d_sk
