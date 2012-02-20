.. -*- coding: utf-8 -*-
.. $Id$

===================
イントロダクション
===================

まずは始めましょう
==================
**PIC code for CANS** （以下 **pCANS** ）はMercurialでバージョン管理されています。MercurialはWindows, linux, Mac OSに対応したバージョン管理ソフトウェアです。Mercurialのインストール・詳しい使い方は、 `本家 <http://mercurial.selenic.com/>`_ もしくは `日本語解説 <http://www.lares.dti.ne.jp/~foozy/fujiguruma/scm/mercurial.html>`_ をご覧ください。

まずは、 **pCANS** のレポジトリを以下のようにダウンロードします。

``hg clone https://bitbucket.org/ymatumot/pic-code-for-cans directory-name``

 **pCANS** は開発途上なため、不定期に更新されます。最新版を反映させるためには、上記で指定した *directory-name* 内で、

``$ hg pull https://bitbucket.org/ymatumot/pic-code-for-cans``

``$ hg update``

とすれば、最新版の差分情報が反映されます。この際、自分で修正を加えたファイルと更新ファイルが重なる場合はマージ（merge）する必要が出てきます。またその結果、衝突（conflict）する可能性もあります。

動作環境
========
- MPI(Message Passing Interface)。2次元以上のコードはMPI並列化バージョンのみです。
- Fortranコンパイラ。ソースコードはFortran 90で書かれています。デフォルトはgfortranです。環境に応じて、Make_incで設定されているMakefileの環境変数$FCを修正してください。
- IDL(Interactive Data Language)。可視化ルーチンはIDLで用意されています。 
- Linuxで動作確認済み。上記ソフト・ライブラリがインストールされていれば、Windows、Mac-OSでもたぶん可だと思います。

環境変数
========
環境変数$PCANS_DIRを **pCANS** がインストールされたディレクトリの絶対パスとし、以下のように設定します（~/pic-code-for-cansにインストールした場合）。

- bashの場合： ``export PCANS_DIR = ~/pic-code-for-cans``
- tcshの場合： ``setenv PCANS_DIR ~/pic-code-for-cans``

可視化について
===============
IDLによる可視化ルーチンが用意されています。$PCANS_DIR/idl内には共通プロシージャ、各課題内にはそれ用のプロシージャがあります。使うためには、環境変数$IDL_STARTUPに$PCANS_DIR/idl/init.proを設定します。

- bashの場合： ``export IDL_STARTUP = $PCANS_DIR/idl/init.pro``
- tcshの場合： ``setenv IDL_STARTUP $PCANS_DIR/idl/init.pro``

$PCANS_DIR/idl/init.pro内には、pathの設定、IDL内の環境の設定等が含まれており、IDL起動時に自動的に設定されます。各自の好みに合わせて修正してください。

IDLについてのさらなる詳細は `こちら <http://www.astro.phys.s.chiba-u.ac.jp/~ymatumot/idl/>`_ を参照ください。

全体の構成
===========
**pCANS** では、1次元及び2次元のコードが用意されています。1次元コードは、シリアル版、MPIによる並列版が用意され、2次元コードはMPI並列版のみとなっています。
