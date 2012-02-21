.. -*- coding: utf-8 -*-
.. $Id$

====================
1次元PICコード
====================


シリアル版
=============
1次元シリアル版コードは以下のような構成になっています。

.. blockdiag::

   diagram {
    "$PCANS_DIR/em1d/" -- "Makefile";
    "$PCANS_DIR/em1d/" -- "Makefile_inc";
    "$PCANS_DIR/em1d/" -- "common/"; 
    "$PCANS_DIR/em1d/" -- "moment/";
    "$PCANS_DIR/em1d/" -- "md_wave/"
    "$PCANS_DIR/em1d/" -- "md_shock/" ;
    "$PCANS_DIR/em1d/" -- "md_whistler/" ;
   }

シリアル版では、物理課題として、「線形波動 (md_wave)」、「衝撃波 (md_shock)」、「電子温度異方性不安定 (md_whistler)」が用意されています（2012年3月現在）。各課題には、

.. blockdiag::

   diagram {
    "md_???" -- "Makefile","main.f90","const.f90","init.f90","dat/","mom/"
   }

が含まれており、それぞれの課題に沿った初期設定のサンプルが置かれています。"dat/"は計算結果の出力先で、電磁場と粒子データが出力されます。"mom/"は"dat/"内にある粒子データを元にモーメント計算をした時の出力先です。

例えば、衝撃波の計算を行うには、

.. code-block:: bash

   $ cd $PCANS_DIR/em1d/md_shock/ 
   $ make
   $ ./a.out

とします。計算終了後にモーメントを計算するには、

.. code-block:: bash

   $ make moment

とします。結果は"mom/"内に出力されます。

パラメタ設定
-------------

MPI並列版
=============
1次元MPI並列版コードは以下のような構成になっています。

.. blockdiag::

   diagram {
    "$PCANS_DIR/em1d/" -- "Makefile";
    "$PCANS_DIR/em1d/" -- "Makefile_inc";
    "$PCANS_DIR/em1d/" -- "common/"; 
    "$PCANS_DIR/em1d/" -- "moment/";
    "$PCANS_DIR/em1d/" -- "md_wave/"
    "$PCANS_DIR/em1d/" -- "md_whistler/" ;
   }

