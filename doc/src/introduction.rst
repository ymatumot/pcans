.. -*- coding: utf-8 -*-
.. $Id$

===================
イントロダクション
===================

まずは始めましょう
==================
**PIC code for CANS** （以下 **pCANS** ）はMercurialでバージョン管理されています。MercurialはWindows, linux, Mac OSに対応したバージョン管理ソフトウェアです。Mercurialのインストール・詳しい使い方は、 `本家 <http://mercurial.selenic.com/>`_ もしくは `日本語解説 <http://www.lares.dti.ne.jp/~foozy/fujiguruma/scm/mercurial.html>`_ をご覧ください。

まずは、 **pCANS** のレポジトリを以下のようにダウンロードします。

.. code-block:: bash

 $ hg clone https://bitbucket.org/ymatumot/pic-code-for-cans directory-name

**pCANS** は開発途上なため、不定期に更新されます。最新版を反映させるためには、上記で指定した *directory-name* 内で、

.. code-block:: bash

 $ hg pull https://bitbucket.org/ymatumot/pic-code-for-cans
 $ hg update

とすれば、最新版の差分情報が反映されます。 **この際、自分で修正を加えたファイルと更新ファイルが重なる場合はマージ（merge）する必要が出てきます。またその結果、衝突（conflict）する可能性もあります。** その場合は出力に従って、text editor等で該当個所を編集してください。

動作環境
========
- Mercurial。レポジトリのダウンロード、最新版のアップデートを行うのに必要です（上述）。
- MPI(Message Passing Interface)。2次元以上のコードはMPI並列化バージョンのみです。
- Fortranコンパイラ。ソースコードはFortran 90で書かれています。デフォルトはgfortranです。環境に応じて、Make_incで設定されているMakefileの環境変数$FCを修正してください。
- IDL(Interactive Data Language)。可視化ルーチンはIDLで用意されています。 
- Linuxで動作確認済み。上記ソフト・ライブラリがインストールされていれば、Windows、Mac-OSでもたぶん可だと思います。

環境変数
========
環境変数$PCANS_DIRを **pCANS** がインストールされたディレクトリの絶対パスとし、以下のように設定します（~/pic-code-for-cansにインストールした場合）。

bashの場合、

.. code-block:: bash

 export PCANS_DIR = ~/pic-code-for-cans

tcshの場合、

.. code-block:: tcsh

 setenv PCANS_DIR ~/pic-code-for-cans

可視化について
===============
IDLによる可視化ルーチンが用意されています。$PCANS_DIR/idl内には共通プロシージャ、各課題内にはそれ用のプロシージャがあります。使うためには、環境変数$IDL_STARTUPに$PCANS_DIR/idl/init.proを設定します。

bashの場合、

.. code-block:: bash

 export IDL_STARTUP = $PCANS_DIR/idl/init.pro

tcshの場合、

.. code-block:: tcsh

 setenv IDL_STARTUP $PCANS_DIR/idl/init.pro

$PCANS_DIR/idl/init.pro内には、pathの設定、IDL内の環境の設定等が含まれており、IDL起動時に自動的に設定されます。各自の好みに合わせて修正してください。

IDLについてのさらなる詳細は `こちら <http://www.astro.phys.s.chiba-u.ac.jp/~ymatumot/idl/>`_ を参照ください。

全体の構成
===========
**pCANS** では、1次元及び2次元のコードが用意されています。1次元コードは、シリアル版、MPIによる並列版が用意され、2次元コードはMPI並列版のみとなっています。$PCANS_DIR内には、以下のようにディレクトリとファイルで構成されています。

.. blockdiag::

   diagram {
    "$PCANS_DIR/" -- "doc/" -- "src/";
    "$PCANS_DIR/" -- "idl/"; 
    "$PCANS_DIR/" -- "em1d/" -- "Makefile","Makefile_inc","common/","moment/","md_???/";
    "$PCANS_DIR/" -- "em1d_mpi/" ;
    "$PCANS_DIR/" -- "em2d_mpi/" ;
   }

"doc/"内には、本マニュアルのソースファイルが含まれています。本ディレクトリに含まれるファイルは開発者向けですので、一般ユーザーは編集する必要はありません。

"idl/"内には、IDLによる可視化ルーチンが含まれています。各課題で使用する可視化のための共通プロシージャが含まれています。

"em1d/"、"em1d_mpi/"、"em2d_mpi/"はそれぞれ、1次元シリアル版、1次元MPI並列版、2次元MPI並列化版コードが含まれます。

各コードのディレクトリ内には、コンパイル用のMakefile、コンパイル時の環境変数を設定したMakefile_incが用意されています。各自の環境によってコンパイラ、コンパイラオプションを指定したい場合は、Makefile_inc内に設定されている、"$FC"と"$FFLAGS"を変更してください。"common/"にはPICコードの共通エンジンが収められています。"moment/"内には、計算結果の粒子データからモーメント計算するためのコードが収められています。"\md_???"は、各物理課題の初期設定等が含まれており、"???"に、物理現象の名前が付けられています。

