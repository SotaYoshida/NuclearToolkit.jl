#!/usr/bin/env python
# coding: utf-8

# # コードの編集環境とGitHub Copilot
# 
# この章では、Pythonを始め様々なプログラミング言語のソースコードを編集するための環境として、VS Codeを紹介する。  
# また、GitHub Copilotという、大規模言語モデルを使ったコード補完ツールを紹介する。
# 
# :::{note}
# 他のノートブックと同様に、Macユーザーは`Ctrl`を適宜`Command`に置き換えて読んでほしい。
# :::
# 
# 免責事項として、本章の記述は2023年7月時点のもの(著者の理解)であり、アップデート等によって内容が不正確となる可能性があることを予めご了承いただきたい。  
# また、一般論として、環境構築に伴うエラーの解決は、使用環境,バージョンや設定などに依存するため、本章の記述をそのまま適用できない場合があること、
# インストール等に関する種々のトラブルについて著者は責任を負わないこともあわせてご了承いただきたい。  
# これも一般論だが、OSに備わっているバックアップ機能(Windows バックアップ, TimeMachine)などで、常時バックアップは取るようにしよう。

# ## Visual Studio Code (VS Code)の導入
# 
# エディタは様々なものがあるが、VS Codeは機能が豊富かつユーザー数が多い(=情報も豊富な)ため、最もオススメしやすいエディタである。  
# とくに拡張機能が便利で、言語ごとに様々な拡張機能が用意されている。著者も紆余曲折ありVS Codeを95%位の割合でメインエディタとして使っている。
# 
# ### インストール方法
# 
# 1. VScodeのインストールは、[公式サイト](https://code.visualstudio.com/)またはOS毎のStore系のソフト(Microsoft Store, Mac App Store)からインストールする。
# 2. インストール後、VScodeを起動すると、日本語化を促されるので、必要であれば指示に従う。
# 3. 初めてPythonのコードをVS Codeで作成/開く際に、拡張機能のインストールを促されるかもしれない。その時は、インストールを許可しよう。
# 4. その他、必要に応じて拡張機能を導入したり、フォントを変更したりすると良い。表示サイズなどはショートカットで拡大(`Ctrl`+`+`)・縮小(`Ctrl`+`-`)できる。
# 
# #### Windows+WSLの導入例
# 
# 著者はWindowsユーザーではないので、情報が古い可能性があるがご容赦頂きたい。
# 
# 1. Storeまたは[Webページ](https://azure.microsoft.com/ja-jp/products/visual-studio-code/)からインストールする
# 
#              
#    <img src="https://github.com/SotaYoshida/Lecture_DataScience/blob/main/notebooks/pic_for_notebook/pic_vscode_1.png?raw=true" width=30%>
# 
# 2.  起動すると、日本語化を提案してくれるので日本語で使いたければ、インストールして再起動をクリック
# 
#    <img src="https://github.com/SotaYoshida/Lecture_DataScience/blob/main/notebooks/pic_for_notebook/pic_vscode_2.png?raw=true" width=30%>
# 
# 3. ターミナルから新しいターミナルを開くと規定のターミナルが開く(Windowsの場合だとpowershell?)
#    ∨記号から、Ubuntu(WSL)を選択すると、Ubuntu(Linux)ターミナルを起動することが出来る
# 
#    <img src="https://github.com/SotaYoshida/Lecture_DataScience/blob/main/notebooks/pic_for_notebook/pic_vscode_3.png?raw=true" width=20%>
# 
#    <img src="https://github.com/SotaYoshida/Lecture_DataScience/blob/main/notebooks/pic_for_notebook/pic_vscode_4.png?raw=true" width=30%>
#    <img src="https://github.com/SotaYoshida/Lecture_DataScience/blob/main/notebooks/pic_for_notebook/pic_vscode_5.png?raw=true" width=30%>
# 
# 
# ### Pythonコードの編集
# 
# 1. 適当なPythonファイルを作成してみよう。[ファイル]→[新しいファイル]で新規ファイルを作成する。  
#    その際**ファイルの種類を入力するか、ファイル名を入力してください**などという表示がでるのでpythonと打ち、Pythonファイルを作成する。  
#    Untitled-1などというファイルが作成されるので、`Ctrl+S`で保存しようとすると、ファイル名を入力するように促される。  
#    その際、適当な名前をつけてみよう。例えば`sample_code.py`などとする。
# 2. ファイルを編集したら、保存(`Ctrl+S`)しよう。編集内容が残っている場合は(標準では)ファイル名の横に丸◯が表示されるので、わかりやすい。
#    Pythonコードだけでなく、様々なファイルを編集したり閲覧したりすることができるので試してみよう。(例: pdfファイルをVS Codeにドラッグ&ドロップしてみよう)
# 3. 拡張子に応じて、なんのプログラミング言語で書かれたソースコードなのかを推定して、色分けしてくれたりする。
# 4. また、VS Codeでは、"フォルダを開く"ことで複数のソースコードを効率的に編集したりすることもできる。
# 
# 
# ### ターミナルの起動
# 
# VS Codeでは、ターミナルをVS Code内で起動することができる。  
# [ターミナル]タブから新しいターミナル(規定のターミナル)を起動することができ、画面の上下や左右などに分割して表示することができる。  
# したがって、ソースコードを編集しながら実行したりといった作業がしやすい。
# 
# Windowsの場合は、WSL(Windows Subsystem for Linux)をインストールしておくと、Linuxのターミナルを使うことができる。  
# MacやLinuxの場合は言うまでもなく、規定のターミナルが開くので、Unix/Linuxコマンドを使えば良い。
# 
# よく使うLinuxコマンドを下記の表にまとめる:
# 
# |コマンド|説明|
# |:--|:--|
# |`ls`|カレントディレクトリのファイル一覧を表示する|
# |`cd`|ディレクトリを移動する|
# |`pwd`|カレントディレクトリのパスを表示する|
# |`mkdir`|ディレクトリを作成する|
# |`rm`|ファイルを削除する(ゴミ箱を経由しないので使用には注意すること)|
# |`cp`|ファイルをコピーする|
# |`mv`|ファイルを移動する|
# |`cat`|ファイルの中身を表示する|
# |`head`|ファイルの先頭を表示する|
# |`tail`|ファイルの末尾を表示する|
# |`grep`|ファイルの中から文字列を検索する|
# |`echo`|文字列を表示する|
# |`chmod`|ファイルのアクセス権を変更する|
# |`sudo`|root権限でコマンドを実行する|
# 
# それぞれのコマンドの使い方やオプションについては網羅的に説明することはしないので、適宜調べて使い方を覚えてほしい。  
# Chat GPTなどに尋ねてみるのも良い。
# 
# ターミナルからPythonコードを実行する際は、
# 
# 1. `python`コマンドでPythonの対話環境を起動して使う
# 2. `python`コマンドに続けてファイル名を入力し、そのファイルを実行する
# 
# の２つの方法がある。最初は1の方法でも良いが、コードで実現したい作業が大きくなるにつれ、
# 作成したソースコードを後者の方法で実行することが多くなるはずだ。
# 
# ```bash
# python sample_code.py
# ```
# 
# ※システムによっては`python3`などとバージョンを明示的に指定しないといけない場合もあるので注意すること。  
# Pythonをコマンドラインで使う際に毎回`python`などと入力するのが面倒なら、
# エイリアスを設定して別ののコマンド(著者は`py`としている)に置き換えてしまうのも良い。
# 
# (Linuxコマンドやシェルについての知識が少しだけ必要になるが)詳しい方法については「python コマンド エイリアス」などで検索してみよう。
# 

# ## 授業資料(.ipynb)をVS codeで実行・編集する
# 
# VS Codeを使うと、ローカル環境で.ipynbファイルを実行・編集することができる。
# 著者も、この授業資料の作成にはVS Codeを使っている。
# 
# 以下で紹介するGitHub Copilotとの相性も良いので、今後プログラミング学習を継続するつもりだ、という方は、
# VS Codeでの.ipynbファイルの実行・編集も是非やってみてほしい。
# 
# ### 授業資料のダウンロード方法
# 
# 本授業資料の.ipynbファイルは、GitHubのレポジトリかBookの各ページからダウンロードすることができる。
# 
# * 方法1: GitHubの[レポジトリ](https://github.com/SotaYoshida/Lecture_DataScience)からダウンロードする.
#   `<>Code`という緑色のボタンをクリックし、`Download ZIP`をクリックすると、全てのファイルをダウンロードできる。
#   展開したフォルダ内の`notebooks`というディレクトリが、本授業資料の.ipynbファイルがある場所である。
#   
#   <img src="https://github.com/SotaYoshida/Lecture_DataScience/blob/main/notebooks/pic_for_notebook/download_repozip.png?raw=true" width=30%>
# 
# * 方法2: Bookの各ページからダウンロードする. ブックの各章の上部にダウンロードボタン(下向き矢印)があるので、そこから`.ipynb`形式を選択し、ダウンロードすれば良い。
#   
# 
# ### ipynb形式のファイルをVS Codeで開く&実行する
# 
# 実行には当然Python環境が必要なので、インストール済みであるとする。
# 
# 1. ターミナルから  
# 
#     ```bash
#     pip3 install jupyter
#     ```
#     などとして、Jupyterをインストールしておく。
# 2. VS Codeを起動し、ダウンロードした`ipynb`ファイルを開く。
# 3. コードセルを実行するには、`Shift+Enter`を押す。または、コードセルの左側にある`▶︎`をクリックする。このとき、複数のPython環境がインストールされている場合は、どのPython環境で実行するかを選択する必要がある。
# 4. コードセルの実行結果が表示される。
# 

# ## GitHub Copilot
# 
# GitHub Copilotは、GitHubがOpenAIと共同で開発したコード補完ツールである。  
# 背後には、GPT-3ベースのOpenAI Codexが動いており、文字通りcopilotとして様々なコードやスニペットを提案してくれる。
# 
# GitHub Copilotは、VS Codeなどの拡張機能としても使用することができる。
# ただし、GitHub Copilotが生成するコードは完璧ではないため、使う際は、**提案されたコードを注意深く確認し、必要に応じて修正を行う必要がある**ことは言うまでもない。
# 
# GitHub Copilotを使うには、まず、GitHub Copilotのベータ版に登録する必要がある。学生は、学生証などを提出することで、無料で登録・使用することができる。上手く使いこなせれば、かなりの生産性向上が期待できるツールであり、私も日々の研究において活用している。
# 
# もちろん、自身の意図を100%反映したコードを生成してくれるわけではないが、コードの骨組みを作ってくれるため大幅に作業が楽になることも多い。ときに「こんな書き方もあるのか」という気づきを与えてくれたりもする。
# 
# 
# ### VS Codeへの導入
# 
# 1. GitHub Copilotに登録する。この際、GitHubアカウントを作成する必要がある。
#    既にGitHubアカウントを持っている場合は、そのアカウントでログインする。  
#    「Github Copilot 学生申請」などで検索すると、学生向けの登録方法が見つかる。
#     認証が済むと、GitHubのアカウントがPROアカウントになり、Copilotを使うことができるようになる。
# 2. 次に、[Copilotのページ](https://github.com/features/copilot/)から、Start my free trialをクリックし、指示に従う。
# 
# 3. VS Codeのアカウント(人型)のアイコンから、GitHubアカウントにログインする。
# 
# 4. VS Codeの拡張機能タブから[GitHub Copilot]などと検索したのち、拡張機能をインストールする(2.でインストールが始まるかも)。その後指示に従い、インストールを完了する。
# 
# ### 試してみよう
# 
# 試しに著者の環境(VS Code)で「適当な二次元データのヒートマップを描画するコードは以下のようになる」などと打つと...  
# 以下のように、自身で入力した際よりも薄い色でコード(場合によっては複数行)が提案される。
# 
# <img src="https://github.com/SotaYoshida/Lecture_DataScience/blob/main/notebooks/pic_for_notebook/pic_copilot_1.jpeg?raw=true" width=70%>
# 
# Copilotの提案を受け入れるには`Tab`キーを押し、棄却する場合には`Esc`を押せば良い。複数行に渡る場合も同様である。
# こうした作業を繰り返して生成されたコードを実行してみよう。
# 

# In[3]:


#適当な二次元データのヒートマップを描画するコードは以下のようになる:

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 2次元データの生成
x = np.random.randn(1000)
y = np.random.randn(1000)

# ヒストグラムの描画
plt.hist2d(x, y, bins=50, cmap='Blues')
plt.show()


# もう少し具体的にカラーマップも指定してみよう。以下のコード例は
# 「適当な二次元データのヒートマップを作成するコード、但しカラーマップはjetを用いる。」と入力したあとに提案されたコードである。

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#データの作成
x = np.arange(0, 10, 0.1)
y = np.arange(0, 10, 0.1)
X, Y = np.meshgrid(x, y)
Z = np.sin(X) + np.cos(Y)

#ヒートマップの作成
plt.pcolor(X, Y, Z, cmap=cm.jet)
plt.colorbar()
plt.show()


# 1回目は`seaborn`ライブラリを使っていたのに対し、２回目では`matplotlib.cm`からカラーマップ`cm`を呼んでいる事がわかる。
# ランダムなデータの作り方も１回目と２回目、あるいはユーザーの実行ごとに変わったりすることに注意しよう。
# 
# 繰り返しになるが、GitHub Copilotの提案が本当に自身の意図したものとなっているかについては、**常にユーザーである我々が注意して確認する必要があること**は肝に命じておこう。

# ### 練習相手・チューターとしてのCopilot/Chat GPT
# 
# この授業を履修する皆さんにはぜひ、GitHub CopilotやChat GPTを、自身の学習のための練習相手や24時間対応のTAやチューターのように活用してもらいたい。  
# 大規模言語モデルを使用したチャットボット(Chat GPT)は、プロンプトと呼ばれる指示文を適切に与えることで、その文脈に沿った文章を生成してくれる(かもしれない)という性質のものである。
# 
# 自分が望むような文章を生成する、あるいは文章生成を"回答"とみなし活用するためには、プロンプトの与え方に工夫が必要になる。  
# 「何がわからないのかも分からない」状態では質問のしようがないのと同様、指示が漠然としているとおそらく有効な回答を引き出す蓋然性は低くなってしまう。
# プログラミング学習に活用するための適切なプロンプトを工夫することは、それ自体が**自身のやりたい作業を言語化する**という作業になっていて、一定の学習効果が期待される。
# 
# 皆さんの主体的な学習にむけて、幾つかChat GPTの使用例を紹介してみよう。  
# 適当な指示文を与えてみる。
# 
# <img src="https://github.com/SotaYoshida/Lecture_DataScience/blob/main/notebooks/pic_for_notebook/pic_copilot_2.jpeg?raw=true" width=50%>
# 
# 
# もちろん、ここでは具体的なファイル構造を示したりファイルを与えているわけではないので、正しく動く保証も全く無いわけだが、概ねもっともらしいコードを生成してくれている。
# 仮に提案されたコードがうまく動かない/使えないとき、うまく行かない理由と関連していそうな情報を与えてやったりすることで、解決策を提示してきたり、補足知識を与えてくれたりするという例も見せよう。
# 
# <img src="https://github.com/SotaYoshida/Lecture_DataScience/blob/main/notebooks/pic_for_notebook/pic_copilot_3.jpeg?raw=true" width=50%>
# <img src="https://github.com/SotaYoshida/Lecture_DataScience/blob/main/notebooks/pic_for_notebook/pic_copilot_4.jpeg?raw=true" width=50%>
# 
# 細かいことを言うと、2.の1つめの画像の提案(`encoding="SHIFT-JIS"`)は、サポートされていないオプションの指定となりエラーになるが、(言語モデルはPythonの文法、とくにPythonの中での文字コードの指定に相当するコードを"理解"しているわけではないため)推論としては真っ当なものになっている。  
# 
# 実は上の画像で挙げた例は、実際にPythonでcsvやエクセルファイルを読み書きしようとする際に、よく遭遇するエラーの例になっていて、
# 私はそのことを知っているため、上のようなある種"誘導的な"プロンプトを与えることで、より有用な回答を得られるよう仕向けている。
# 2.の2つめの画像で、より適切なオプションを指定するためにもう少しヒントをあげているという訳だ。
# 
# 
# "自身の理解やこれまでの思考・試行を開示した上で質問をする"ことが、他人に質問する際に有効な回答を引き出すための近道であることと同じで、
# 言語モデルを対話的に活用する上では、こうしたコツを抑えておくと良いかもしれない。
# 
# また、大規模言語モデルの性質・特性やその背景にある基礎について理解を深めることも非常に重要な姿勢だと感じる。  
# 理研AIPシンポジウムでの岡崎直観先生の特別講演[Youtubeへのリンク](https://youtu.be/PUuk4Cv-ycg)が非常に参考になる。ぜひ一度視聴してみてほしい。
