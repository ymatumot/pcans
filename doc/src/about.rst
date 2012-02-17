.. -*- coding: utf-8 -*-
.. $Id$

====================
reSTの書き方について
====================

このドキュメントは `sphinx`_ ( `日本語版ドキュメント`_ )で生成されています。sphinxは `reStructuredText`_ (reST)という形式(+いくつかの拡張)で書かれたテキストファイルをhtmlやpdfへ変換するツールです。 `数式の表示`_ や `ソースコードのハイライト表示`_ も可能です。

.. _sphinx: http://sphinx.pocoo.org/
.. _日本語版ドキュメント: http://sphinx-users.jp/doc.html
.. _reStructuredText: http://docutils.sourceforge.net/rst.html

数式の表示
==========

以下のように簡単にLaTeX形式で数式を書くことができます。
まずはディスプレイスタイルでMaxwell方程式を表示します。

.. math::

   \frac{1}{c}\frac{\partial \mathbf{B}}{\partial t} &= -
   \mathbf{\nabla} \times \mathbf{E} \\
   \frac{1}{c}\frac{\partial \mathbf{E}}{\partial t} &=
   \mathbf{\nabla} \times \mathbf{B} - \frac{4 \pi}{c} \mathbf{J}


文中に :math:`\lim_{x \rightarrow 0} \frac{\sin x}{x}` のようにインラインで埋め込むことも可能です。LaTeXでいうところの ``$数式$`` です。

ソースコードのハイライト表示
============================

以下のように自動的にソースコードのハイライト表示ができます::

   program test
   implicit none
   integer*4 i

   write(*,*) 'Hello World !'

   do i = 1, 10
      write(*,*) i
   end do

   end program test

他の言語でも可能です。以下はC言語の例です。

.. code-block:: c

  #include <stdio.h>

  int main()
  {
   printf("Hello World !\n");
   return 0;
  }

また別のファイルをincludeして表示することもできます。
