.. -*- coding: utf-8 -*-
.. $Id$

===================
イントロダクション
===================

インストール
=============
**PIC code for CANS** はMercurialでバージョン管理されています。MercurialはWindows, linux, Mac OSに対応したバージョン管理ソフトウェアです。Mercurialのインストール・詳しい使い方は、 `本家 <http://mercurial.selenic.com/>`_ もしくは `日本語解説 <http://www.lares.dti.ne.jp/~foozy/fujiguruma/scm/mercurial.html>`_ をご覧ください。

まずは、PIC code for CANSのレポジトリを以下のようにダウンロードします。

``hg clone https://bitbucket.org/ymatumot/pic-code-for-cans directory-name``

PIC code for CANSは開発途上なため、不定期に更新されます。最新版を反映させるためには、上記で指定したdirectory-name内で、

``$ hg pull https://bitbucket.org/ymatumot/pic-code-for-cans``
``$ hg updatehg pull``

とすれば、最新版の差分情報が反映されます。この際、自分で修正を加えたファイルと更新ファイルが重なる場合はマージ（merge）する必要が出てきます。またその結果、衝突（conflict）する可能性もあります。mergeする場合は、Commitsタブで最新版の変更点を確認してください。 

