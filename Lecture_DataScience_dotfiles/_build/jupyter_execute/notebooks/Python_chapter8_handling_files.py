#!/usr/bin/env python
# coding: utf-8

# <a href="https://colab.research.google.com/github/SotaYoshida/Lecture_DataScience/blob/main/notebooks/Python_chapter8_handling_files.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# # ファイル・文字列操作
# 
# [この章の目的]
# text,csvやxlsx形式のデータをプログラムでサクッと扱えるようになる。

# この章では、テキストファイルやcsvファイル(excelファイルはおまけ$\clubsuit$)をPythonで操作する簡単な方法を学習する。  
# 
# これまでの章では、データはリストとして既に与えられた状態から解析を行ったが、実際にデータを扱う際は  
# 既に誰かが作成した何らかのファイルをプログラムで読み込んで操作する場合も多い。  
# この章の内容は、データ解析というよりは、Pythonでデータ解析をするための下準備に相当する。  
# 
# 愚直にコードを書いている事もあり少々泥臭い部分が多いが、この章のような操作のエッセンスを抑えておけば  
# 普通にやると膨大な時間がかかる様々な処理を高速化・自動化することができるので、頑張って学習しよう。
# 

# ## 授業で使うファイルの準備
# 
# 予め以下のリンクをクリックして、ファイルをダウンロードし、  
# ご自身のGoogle Driveにアップロードしておいてください。
# 
# * [test.txt](https://drive.google.com/file/d/1U2uvrN18713ylN4OQiI2fsfX5gudL45w/view?usp=sharing) (テキストファイル)
# 
# * [python_handling_test.csv](https://drive.google.com/file/d/1bYJNWdtujcQWfSBAa1UeXi2ZzJRJktil/view?usp=sharing) (csv, カンマ区切りのテキストファイル)
# 
# * [kakei.xlsx](https://drive.google.com/file/d/1gJMVHivmP7R9Qf4LdqRhdPVc3x0IzD8v/view?usp=sharing) (エクセルファイル)
# 
# 本章では、ファイルの場所を指定する**パス**という概念がたびたび登場する。
# 以下のコードをそのまま使いたいという方は、マイドライブ直下に`AdDS`というフォルダを作り、さらにその中に```chapter8_data```というフォルダを作成し、ファイルをいれてください。
# 
# パスについては後の節で詳しく説明します(今気になる方は末尾をチェック)。

# ## テキストファイルの操作

# 膨大な行のテキストファイルに対して、人間が手で操作をするというのは時として非現実的です。  
# 
# 誤変換を置換するくらいなら、どのテキスト/メモ帳アプリやwordでもできますが、
# 全行(数千とか数万)に対して、決まった操作が必要な場合、プログラムにしてしまったほうが遥かに便利です。
# 
# 以下ではGoogle Driveのマイドライブの```AdDS```の下に作った```chapter8_data```というフォルダにファイルを保存したと仮定して話を進めますので、**適宜皆さんの場合に置き換えて使用してください**
# 

# まずはgoogle driveに保存した```test.txt```という名前のファイルを読み込んでみましょう。  
# 既に何回かやったようにgoogle driveをマウントします。
# 
# 
# 

# In[ ]:


from google.colab import drive
drive.mount('/content/drive')


# **注意** 以後のコードは、google driveの中にあるファイルを読み書きしたりといった操作を行うため  
# 上でGoogle Driveをマウントした状態でなければ実行しても多くがエラーとなる。
# 
# ---
# 
# Google Driveのマウントができたら、先程のファイルがあるかlsコマンドで確かめてみよう。
# 

# In[ ]:


get_ipython().system('ls /content/drive/MyDrive/AdDS/chapter8_data/*')


# *はワイルドカード記号で、対象を任意とする命令に相当します。  
# *.拡張子 などとして使うのも便利です。

# ファイルが見つからない場合は
# > No such file or directory
# 
# などと表示される。

# > $\clubsuit$ 上のファイルをGoogle Driveに保存する方法としては  
# 一度ローカルに保存したファイルをブラウザなどからアップロードする方法はもちろん  
# 全てを(Linux)コマンドで行うこともできる。
# ```
# !git clone https://github.com/SotaYoshida/Lecture_DataScience
# !mkdir /content/drive/MyDrive/AdDS/
# !mv Lecture_DataScience/Chapter8_data /content/drive/MyDrive/AdDS/chapter8_data
# !ls /content/drive/MyDrive/AdDS/chapter8_data
# ```
# 1つめの行ではまず授業資料のGitHubレポジトリをColab環境で間借りしているgoogleのサーバー上にクローン(≒コピー)する。2行目でマイドライブの下にAdDSというフォルダの作成を試み、3行目でダウンロードしてきたレポジトリにある`Chapter8_data`をさっき作ったAdDSというフォルダの中に別名(先頭が小文字になっている)で移動する。
# 最後に、どんなファイルがあるかをlsコマンドで確認している。  
#   重複する作業も多いのでこれらのコードはコードセルには書かなかったが、うまくアップロードできなかった場合やコマンドラインによるファイルの移動などをやってみたければ上の4行のコードをコードセルに貼って試してみよう。

# In[ ]:


get_ipython().system('ls /content/drive/MyDrive/AdDS/chapter8_data/*txt ')


# とするとマイドライブ/AdDS/chapter8_data/にある`.txt`形式のファイル一覧を表示させることができる。  
# `test.txt`が見つかったでしょうか？(DriveにアップロードしてColabから読み込みできるまでに少し時間がかかる場合がある)
# 
# 
# では次に、このファイルに書かれているテキストを取得してみよう。  
# 方法は幾つかあるが、最も標準的なものとして、ファイルを開いてテキストを取得する方法を試してみよう。

# ### テキストファイルを開いて内容を読み出す

# In[ ]:


filename = "/content/drive/My Drive/AdDS/chapter8_data/test.txt" 
inp = open(filename,"r")
lines = inp.readlines()
inp.close()


# 1行目でファイル名(正確にはファイルのパス)を指定し```filename```という変数にした。  
# 
# 2行目では、指定したパスにあるファイルを開いている。  
# 今はファイルに書き込むのではなく、既にあるファイルを開いて読み込むので`"r"`というオプションを指定している。  
# 他には`"w"`(書き出し,上書き), `"a"`(書き出し,追記)などがあり、新しく上書きでファイルを作成したい場合は`"w"`,すでにあるファイルの内容は消さずに追記したい場合は`"a"`を指定して使う。
# 
# 3行目では、`inp`(ファイルを`open`して得たオブジェクト)に対して```readlines```という操作を適用している。
# これは、ファイルに書かれているテキストを(可能なら)全行に渡って読み込みメモリにストアする関数になっている。  

# In[ ]:


print(lines)


# とすると、全ての行が読み込まれ、変数```lines```に格納されていることがわかる。ここで```\n```は改行記号を意味する。
# 
# ループを回して一行ずつ表示させると

# In[ ]:


for line in lines:
    print(line)


# といった感じ(行ごとにスペースが生じている理由については後で説明します)。

# 必要な行番号が分かっている場合は、`nlines = lines[2:10]`などとして、要らないところは捨てても良い。(リストのスライスについては2章を参照)
# 
# 次に、もう少し具体的なテキスト操作をしてみよう。   
# まず、上の1行ずつ表示するコードでは、改行コードを明示的に含む文字列を一行ずつ表示したため、改めて`print`すると余分なスペースが空いてしまう。  
# (`print`関数はデフォルトで末尾に改行```\n```を挿入するのでファイルにある改行記号とあわせて2回改行してしまう→[参考リンク](https://docs.python.org/ja/3/library/functions.html#print))

# ### strip関数
# 
# たとえば```strip()```関数を使うと、文字列に含まれる空白、タブや改行コードを消去することができる。

# In[ ]:


a = "test character\t"
b = "test2 \n"
print("a", a, "←タブが隠れている")
print("b", b, "←改行される")
### strip関数をもちいて...
print("a.strip()", a.strip(),"b.strip()",b.strip())


# 先程のforループでstrip関数を適用してやると...

# In[ ]:


for line in lines:
    print(line.strip())


# 文字列の右側に空白や改行コードが入っていることが明確な場合は  
# `strip`の代わりに`rstrip`を使ってもOK(`rstrip`のrはrightの意味)。

# In[ ]:


for line in lines:
    print(line.rstrip())


# 
# ファイルによってはインデントをするために左側にタブ```\t```が含まれる場合もあります(PythonのコードをテキストとしてPythonから読むときなどがこれに該当)。そのような場合に左側にある空白やタブのみを取り除きたければ```lstrip()```を使って取り除くことができる。
# 
# もちろんPythonではインデントが文法なので、インデントを一律で消す、といった操作は必要ないが、特定の状況では、`lstrip`も使えると便利だ。

# 上のファイルの文字列で`#`記号を含む行以降だけが必要な場合はどうすればいいでしょうか？
# 
# 最も単純(？)な実装は

# In[ ]:


hit = 0 #
for line in lines:
    if "###" in line:
        hit += 1 
        continue
    if hit == 0 :
        continue #hitが0の状態では何もしない
    print(line.rstrip())


# といった愚直な例が考えられる。
# つまり、`#`を含む行に到達するまでの行は無視して、それ以降の行だけをリストに格納するというもの。
# もちろん`#`を含む行が複数あるようなケースでは、自分が実現したい操作にあわせてコードを書き換える必要がある。
# 
# 以下では、`###data`までの行が必要ないので、必要なところまでを別のリストに格納してしまおう。
# 

# In[ ]:


hit = 0 #
nlines = []
for line in lines:
    if "###" in line:
        hit += 1 
        continue
    if hit == 0 :
        continue #hitが0の状態では何もしない
    nlines += [line]
print(nlines)


# ### split関数

# また、1,2,3,4,5,6といったコンマやスペースで区切られたものをリストに格納したい場合には、```split```関数が便利。`split`関数は引数に何も指定しなければ、スペースや改行もしくはタブごとに文字列を区切ったリストを返す。

# In[ ]:


sample_text = "This is a\nsample\ttext."
sample_text.split()


# In[ ]:


for line in nlines:
    print(line.split())


# カンマがあるときはカンマで分割する、という約束を表現したければ

# In[ ]:


for line in nlines:
    if "," in line :
        print(line.rstrip().split(","))
    else :
        print(line.rstrip().split())


# などとすれば良い。これを利用すれば、空のリストにファイルから読んだ要素を詰めていくといった操作も実現できる。
# 

# In[ ]:


# 数字とプロフィールの空リストを作り、そこに読み込んだものを詰めていく
# その際に、数字のリストとプロフィールのリストを分けたいとする
nums = [] 
profs = [] 

for line in nlines:
    if "," in line :
        nums += [ line.rstrip().split(",") ]
    else :
        profs += [ line.rstrip().split()]
print("nums", nums)
print("profs", profs)


# 上の`nums`の様に、予め全ての要素が整数だと分かっていて整数に対する演算(四則演算など)を後でするのなら、`str`(文字列)型ではなく`int`型にしておくほうが良いこともあるだろう。

# In[ ]:


##リスト内包表記を使った実装
nums = []
for line in nlines:
    if "," in line : 
        tl =  line.rstrip().split(",")
        nums += [ [ int(tmp) for tmp in tl] ]
print("方法1:", nums)

## map関数(後述)を使った実装
nums = []
for line in nlines:
    if "," in line : 
        nums += [ list(map(int, line.rstrip().split(",") )) ]
print("方法2:", nums)


# ### replace関数

# `replace`関数で文字の置換が可能です

# In[ ]:


text = "abcdあいうえお"
text = text.replace("abcd", "1234")
print("置換や→",text)
print("除去にも→", text.replace("4", ""))


# ### $\clubsuit$ map関数
# 
# `map`関数は`map(操作,対象)`という風に使って、対象の各要素に対して一括で操作を適用することができます。  
# 今の場合、`['1', ' 2', ' 3', ' 4', ' 5', ' 6']`などの文字列のリストに対して、  
# 整数型に変換する```int```関数を作用させるという操作を一度に行います。
# 
# >注: `map`関数の返り値はmap objectと呼ばれるものなので、  
# 元のようなリストの形で使いたい場合は```list()```を使ってリストに変換するステップが必要です。
# 
# 世の中には、アンケート結果や産業データなどがcsv(カンマ区切りのテキスト)ファイルで公開されている場合が多いですが、  
# その場合は**上で説明したような手順でリストなどに格納すれば今まで行ったような解析やグラフ描画が実行できる**といったわけです。  
# (もちろんcsvを読むのに便利なライブラリもありますが  
# いろんな形式のファイルをプログラムで読み込む場合には  
# 上のような基本的な操作を組み合わせることも必要になります。)

# ### テキストファイルの書き出し

# 次に、テキストファイルを書き込んで保存してみます。  
# 上の文字列で、敬称を"さん"から"様"に置換したテキストを作成して、それを別ファイルとして保存してみましょう。
# 

# In[ ]:


filename = "/content/drive/My Drive/AdDS/chapter8_data/test_replace.txt" 
oup = open(filename,"w")  ## oup は"output"の気持ち...
for line in lines:
    print(line.rstrip().replace("さん","様"), file=oup) # file=[openしたファイル]にすることで、printする先をファイルに指定できます。
oup.close() #ファイルはきちんと閉じる.


# Google Driveで、作成されたファイルをチェックしてみましょう。
# 
# なお、filenameに元ファイルと同じものを指定すると```open(filename,"w")```を実行した時点で  
# ファイルが上書きされて空ファイルになるので注意しましょう。

# 今の例ではもちろん、手で置き換えたりするほうが遥かに速いですがこうしたPythonによるファイル操作を覚えておくと
# 
# * ファイル自体が大量にあり、同じ操作を繰り返す場合
# * 単一のテキストファイルに大量の行に渡って情報がある場合
# 
# など、手作業が非現実的な様々な状況でも、楽ちんに作業を終わらせることができる(かもしれません)。
# 
# 上の内容や、これまでに学習したループ処理を駆使すると、  
# 数万人のデータが1行ずつ記載されたテキストファイルから条件にヒットする人の  
# 情報だけを抽出して小さなファイルにまとめるといったことも可能です。
# 
# **プログラミングを用いたファイル操作をする発想があるかどうか**がきっとこの先  
# 皆さんの生き方や働き方に大きな影響を与えると私は考えています。
# 
# > **文字コードに関連した余談**  
# Windows環境で作成されたテキストファイルを扱う際は読み込みで、文字コードによるエラーが出るかもしれない。最近ではメモ帳でもUTF-8(世界標準)を採用しているよう(→[MicrosoftのWindows blogの記事](https://blogs.windows.com/japan/2020/02/20/about-windows-and-japanese-text/))だが、古いテキストファイルだとShift-JISになっているかも。そういうときは、```open(file, "r", encoding = "shift_jis")```など、ファイルを開くときにencodingを明示的に指定する必要がある。明示的にUTF-8で保存したいときは```open(file, "w", encoding = "utf-8")```などとする。  
# 参考: [公式ドキュメント](https://docs.python.org/ja/3/howto/unicode.html#reading-and-writing-unicode-data)  
# ここまで勉強してきた皆さんには「そんなの、パソコンに存在するShift-JISで書かれたテキストファイルを全てUTF-8に変換するPythonスクリプト書けばいいんじゃね？」という発想があることを期待しています。
# 

# ## csv,エクセルファイルの操作

# ### アンケート分析
# 
# 冒頭の二番目のファイル[python_handling_test.csv](https://drive.google.com/file/d/1bYJNWdtujcQWfSBAa1UeXi2ZzJRJktil/view?usp=sharing)はあるアンケート結果をまとめたファイルになっています。
# 
# これは、Google フォームで作成したアンケーで、国数英社理(中学の５科目)に対する得意/苦手意識の調査を想定した疑似アンケートです。
# 
# このようなアンケート調査は事務作業や卒業研究などで頻繁に見られ、会社や大学など所属コミュニティで何らかの意思決定に用いられることも多いことでしょう。こうしたアンケート分析を行っていると、
# *   各回答項目同士の関係が知りたい
# *   明示的な項目以外の情報も抽出したい
# 
# といった要望が出てきます。今の場合でいうと、
# * 各科目ごとの得意・苦手意識の相関を調べたい
# * 夜中(あるいは日中)にアンケートを回答した夜型(昼型)の人に見られる特徴がなにかないか？
# 
# といったイメージです。そんなとき、
# 
# > 国語が得意(どちらかというと得意)と回答した方に質問です。  
# 英語についてはどうでしょうか？
# 
# などと新たに設問を増やしてアンケートをやり直すというのは得策では有りません。  
# すでに得られた情報からさらなる情報を引き出すことを考えてみましょう。  
# まずは、csvファイルに記載された情報を整理してプログラムで扱いやすくすることを考えます。
# 
# > 余談: このcsvファイルをExcelで開こうとするとお使いの環境によって文字化けを起こすかと思います。これはgoogleフォームで作成されたcsvファイルの文字コードが世界標準のutf-8を使用しているのに対し、ExcelがShift-JISという時代遅れな文字コードでcsvファイルを開こうとするためです。Googleのスプレッドシートや、Mac標準のNumbersで開くと文字化けしません。
# 
# > 2000件の回答は、もちろん私が手作業で入力したわけでも誰かに協力してもらったわけでもなく、一定のルール(傾向)を勝手に設定した上でランダムに回答を作成しフォームから自動回答するPythonスクリプトを書きました。  
# 時間に余裕があれば、こうしたWeb操作を自動化する方法も授業で扱います。 c.f. ブラウザ操作, Webスクレイピング

# In[ ]:


filename = "/content/drive/My Drive/AdDS/chapter8_data/python_handling_test.csv" #読み込むファイルのパスの指定


# とりあえずファイルの中身を数行表示してみる。
# 
# csvファイル(コンマ区切りのテキスト)なので、テキストファイルと同じ方法をとる(他の方法ももちろんある)

# In[ ]:


inp=open(filename,"r")
csv_lines=inp.readlines() 
inp.close()
print("行数は",len(csv_lines))
for i in range(5):
    print(csv_lines[i].rstrip())


# ちなみに...```pandas```ライブラリを使うとcsvをサクッと読み込むことができる

# In[ ]:


import pandas as pd 
df = pd.read_csv(filename)
print(df)


# さて、```csv_lines```に格納したデータをもう少し扱いやすいように変更しよう。  
# 最初の0行目はどういうデータが入っているか(データの項目)を表している。  
# 1-2000行目には2000人分の回答が詰まっている。  
# 
# これによると、  
# > 0列目: 回答した時刻  
# > 1列目: 性別  
# > 2列目: 国語  
# > 3列目: 数学  
# > 4列目: 英語  
# > 5列目: 社会  
# > 6列目: 理科  
# 
# らしい。いろいろなデータの整理方法があると思うがここでは、
# * 処理A 0列目の時刻を24時間表記にして表示する  
# * 処理B 2-6列目の各科目の得意・苦手意識を、文字列を除去して数値[-2,-1,0,1,2]として扱う
# 
# をまずやってみよう。
# 
# 
# 
# 
# 

# In[ ]:


#処理Aのための関数
#input_strが、"年月日 時刻(h:m:s) 午前/午後 GMT+9" という文字列である、という文字列の[構造]を知った上での実装になっていることに注意
def make_time_24h(input_str):        
    time  = input_str.split()[1]
    AMPM = input_str.split()[2]
    hms = time.split(":")
    h = int(hms[0])
    if AMPM == "午前":
        output_str = time 
    else :
        if h != 12:
            output_str = str(h +12)+":"+hms[1]+":"+hms[2]
        else:
            output_str = str(h)+":"+hms[1]+":"+hms[2] # 12時xx分だけは別の取り扱いが必要
    return output_str

nlines=[] #整理したものをリストとしてまとめるための空のリスト
for nth,line in enumerate(csv_lines[1:]): 
    nline = line.rstrip().replace('"','').split(",") # 改行文字の除去、ダブルクォーテーションの除去, カンマで分割    
    # この時点でnlineは0:時刻 1:性別, ...のリストとなっているはず print()でcheckしてみよう
    # 処理A)
    time = make_time_24h(nline[0])
    #print("nline[0]", nline[0], "time", time)
    M_or_F = nline[1] #性別

    #　処理B)
    points = [ int(nline[k].split()[0]) for k in range(2,7)] #各科目の値だけのリスト(points)を作成
    # 上記をmap関数にしてみよう。

    nline = [time, M_or_F]+points  #リストを連結(時刻,性別と各科目の値を同じ階層で結合)して、nlineという名前で上書き
    nlines += [ nline ]

    # うまく編集できたか400行おきほどでprintしてチェックしてみる
    if nth % 400 == 0 :
        print("編集前", line.rstrip())
        print("編集後", nline)
        print("")


# 最後に、各項目の得点を適当なリスト(あるいはnp.array)に整形しておけば、種々の分析を行うことができます。
# 
# 

# In[ ]:


import numpy as np
points = [ [] for i in range(5)]
for tmp in nlines:
    for i in range(5):
        points[i]+=[tmp[2+i]]
print("points", np.array(points))
print("各科目の平均スコア:", [np.mean(points[i]) for i in range(5)])


# 相関分析は以降の章で扱うので具体例は省略します。
# 

# ## $\clubsuit$ 複雑なエクセルファイルの操作
# 
# ```kakei.xlsx```はエクセルファイルで以降では、2020年度前期のデータサイエンス入門(一部学科を除く)の  
# 相関分析で使用されたエクセルファイル、[kakei.xlsx](https://drive.google.com/file/d/1gJMVHivmP7R9Qf4LdqRhdPVc3x0IzD8v/view?usp=sharing)を使用します。  
# 
# 
# 以下では、上と同じディレクトリに`kakei.xlsx`を置いたと仮定して  
# 処理を行いますので、適宜ご自身の環境にパスを置き換えてください。
# 
# ※もともとはxlrdというライブラリを使って実装していましたが.xlsx形式をサポートしなくなるとのことで、pandasライブラリを用いた実装に変更しました。

# In[ ]:


filename = "/content/drive/My Drive/AdDS/chapter8_data/kakei.xlsx" #読み込むファイルのパスの指定


# まずはxlsxファイルをPythonで読み込んで、どんな"シート"があるのかを確認してみましょう。

# In[ ]:


import pandas as pd
input_file = pd.ExcelFile(filename)
sheet_names = input_file.sheet_names
print("pandas: シート名",sheet_names)


# たくさんシートがあることが分かります。 
# 

# Sheet1の中身をのぞいてみましょう。まずは行と列の数を取得してみます。

# In[ ]:


Sheet1 = pd.read_excel(filename, sheet_name="Sheet1")
print("行,列の数", Sheet1.shape)


# 0-5番目の行にはどんな値がセルに入っているのかな...と思ったら

# In[ ]:


for i in range(5):
    print( list(Sheet1.iloc[i]) )


# などとする。このように、扱いたいファイルの"構造"を知ることが  
# やりたい操作を系統的に実行するための第一歩です。  
# このエクセルを実際に開くとSheet1からSheet12までが複数都市の家計調査のデータで  
# S1からS12までが気候データになっていて  
# 1-12までの数字が2017年の1月から12月までに対応していることが分かります。
# 
# 実際のデータを触っていると「2006年までとそれ以降とでデータファイル(.xlsx)の"構造"が違う」  
# といったことも出てきます。  
# 最初は特定のものに合わせたコードを作り、徐々に"汎用性の高い"コードにしていくのがよいでしょう。
# 
# このエクセルを使って実際に作業をするには、[細かいライブラリの使い方]などを説明することになるため   
# 授業ではやらず、以下の"おまけ"にいれておきます。この作業や実践DSに限らず
# * 自分がやりたい操作をきちんと言語化する
# * 公式ドキュメントやWebから情報を探す
# * とにかく試してみる
# 
# という意識が重要です。
# 

# ### $\clubsuit$$\clubsuit$おまけ
# 
# 以下のコードは、プログラミングの"ありがたみ"を感じてもらうためのお試し用です。   
# (昔書いたかなり読みにくいコードなのであまり真剣に読まないでください.)
# 
# **大量の画像ファイルをドライブに生成するので、以下を読んだ上で実行してください**
# 
# 以下のコードたちを何もいじらずに実行すると、  
# 全都市の月別平均気温と全品目の世帯平均支出のうち、  
# 相関係数の絶対値が0.9以上のもの(291通り)をプロットして画像として保存します。  
# ```pthre```の値を小さくすると、生成される画像の数がとんでもなく増えるのでやらないでください。
# 
# (0.9 →  291通り, 0.8 → 1234通り, 0.7 → 2871通り,  
#  0.6 → 5233通り, 0.5 → 8375通り, 0.0 → 32876通り)
# 
# Google Colab上で実行して291枚の画像が生成されるまでに80~150秒程度かかるようです。
# 
# この時間未満でエクセルで操作をして同様の処理を完了出来るという方は...おそらく地球上にいないでしょう(要出典)
# 
# 

# In[ ]:


# 画像がいっぱい生成されると面倒なので画像を保存するフォルダを作成しておく
get_ipython().system('mkdir /content/drive/MyDrive/AdDS/chapter8_data/kakei_cor_pic ')


# In[ ]:


get_ipython().system('pip install japanize_matplotlib ')


# In[ ]:


import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import japanize_matplotlib
import time

class ebook:
    def __init__(self,inpf):
        self.input_file = pd.ExcelFile(filename)
        sheet_names = input_file.sheet_names
        self.sname = sheet_names
        self.ns = len(sheet_names)
        print("pandas: シート名",sheet_names)
        print("self.ns", self.ns)

        s_kikou=[]; s_kakei=[]
        for i, sheetname in enumerate(self.sname):
            if "Sheet" in sheetname :
                s_kakei += [ i ]
            elif "S" in sheetname :
                s_kikou += [ i ]
        self.s_kakei,self.s_kikou = s_kakei,s_kikou
    def indices(self):
        return self.s_kakei, self.s_kikou
    def readkakei(self,ikakei) :
        ws = self.input_file.parse(sheet_name=self.sname[ikakei])
        nr = ws.shape[0]
        premode = True
        items = []
        for ii in range(nr): 
            trow = list(ws.iloc[ii])
            hit = 0
            if premode == True:
                for jj,tmp in enumerate(trow):
                    if type(tmp) is str:
                        if "市" in tmp:
                            hit += 1
                if hit > 5:
                    premode=False
                    i_kakei=[];p_kakei=[]
                    for jj,tmp in enumerate(trow):
                        if type(tmp) is str:
                            if "市" in tmp:
                                i_kakei += [jj]
                                p_kakei +=[ tmp ] 
                    v_kakei = [ ]
            else:                    
                if ii >= 22:
                    if type(trow[8]) is str and trow[8] != "":
                        v_kakei += [ [trow[jj+1] for jj in i_kakei] ]
                        items += [trow[8]]                         
        return i_kakei, p_kakei, v_kakei,items
    def readkikou(self,ikikou):
        ws = self.input_file.parse(sheet_name=self.sname[ikikou], header=None)
        nr = ws.shape[0]
        quantities = [];v_kikou=[]
        premode=True
        for ii in range(nr): 
            trow = list(ws.iloc[ii])
            if premode :
                if any(["市" in str(tmp) for tmp in trow]):
                    Tplaces = trow[1:]
                    premode=False
            else:
                quantities += [ trow[0] ]
                v_kikou += [ trow[1:] ]
        return Tplaces, v_kikou,quantities

def seasoncolor(month):
    if month <= 2 or month ==12:
        return "blue"
    elif 3 <= month <=5:
        return "green"
    elif 6 <= month <=8:
        return "red"
    elif 9<= month <=11:
        return "orange"
    return tcol

def plot_cor(x,y,item,quantity,place,corrcoef):    
    fig = plt.figure(figsize=(4,4))
    ax = fig.add_subplot(1,1,1)
    ax.set_facecolor("#e0e0e0")
    ax.set_title(place+"   r="+str("%8.2f" % corrcoef).strip())
    ax.set_xlabel(item);ax.set_ylabel(quantity)
    ax.grid(True,axis="both",color="w", linestyle="dotted", linewidth=0.8)
    for i in range(len(x)):
        tcol=seasoncolor(i+1)
        ax.scatter(x[i],y[i],marker="o",s=5,color=tcol,zorder=20000,alpha=0.7)
        ax.text(x[i],y[i],str(i+1)+"月",color="k",fontsize=8)
    plt.savefig(oupdir + "corr_"+item+"vs"+quantity+"_at_"+place+".png",dpi=300) 
    plt.close()

def calcor(date,places,items, Vs_kakei,Tplaces,quantities,Vs_kikou):
    hit = 0; num_pic=0
    Vs = [] 
    for j_K,place in enumerate(places):
        for j_T, Tplace in enumerate(Tplaces):
            if place != Tplace :
                continue
            for ik,item in enumerate(items):
                kvalue = np.array([ Vs_kakei[i][ik][j_K] for i in range(len(Vs_kakei))])
                quantity=quantities[iT]
                Tvalue = np.array([ Vs_kikou[i][iT][j_T] for i in range(len(Vs_kikou))])
                if all(Tvalue) == 0.0: ## missing value in climate data
                    continue
                if printlog:
                    print("@", place," ",item,kvalue," VS ",quantity, ",",Tvalue)
                corrcoef=np.corrcoef(kvalue,Tvalue)[0][1]
                Vs += [ [ corrcoef, item, quantity, place] ]
                if abs(corrcoef) > pthre:
                    hit += 1
                    if pltmode==True:
                        plot_cor(kvalue,Tvalue,item,quantity,place,corrcoef)                       
                        num_pic += 1
    print("hit:",hit, " number of picture", num_pic)

if __name__ == "__main__":
    ti=time.time()
    T=True; F=False

    inpf = "/content/drive/My Drive/AdDS/chapter8_data/kakei.xlsx"
    oupdir = "/content/drive/My Drive/AdDS/chapter8_data/kakei_cor_pic/" #適宜置き換える
    iT = 6  # iT=6: 日平均気温
    printlog= F #条件にhitした都市の品目と気候データを逐次printするかどうか. (Fを推奨)
    pthre= 0.90 ## corrplotを描く相関係数のthreshold 
    pltmode = T ## T:plotする F:計算のみ ** 画像をいちいちplotして保存する必要がない場合Fを推奨
    year="2017" 

    wb=ebook(inpf)
    s_kakei,s_kikou=wb.indices()   
    Vs_kakei=[]; Vs_kikou=[];dates=[]
    for i,ind_kakei in enumerate(s_kakei):
        i_places,places, v_kakei,items = wb.readkakei(ind_kakei)
        Tplaces, v_kikou, quantities  = wb.readkikou(s_kikou[i])
        if i+1 < 10:
            date=year+"0"+str(i+1)
        else:
            date=year+str(i+1)
        dates += [date]
        Vs_kakei += [ v_kakei ]
        Vs_kikou += [ v_kikou ]
    calcor(dates,places,items,Vs_kakei,Tplaces,quantities,Vs_kikou)    

    tf=time.time()
    print("Elapced time[sec]:", tf-ti)


# ## 余談: 電子ファイルのフォーマット
# 
# プログラムでデータを機械的に読み出して活用することで、人間が到底出来ないような作業効率を実現することができる場合も多い。
# そんな光の側面ばかりなら良いが、実際にはそう上手くは行かないことも多い。
# 
# 業務のデジタル化・デジタルトランスフォーメーションなどといった標語とは裏腹に、世の中にあふれるcsv,スプレッドシートなどは、
# csvと謳っておいて、実際にはカンマ区切りではなくタブ区切りであったり、機械判読を全く想定していないデータの書き方・並べ方となっているものが多く、プログラムを書ける人にとっては苦痛な状況も多い。  
# 
# 総務省統計局は令和2年2月に、政府統計(e-Stat)に関して[統計表における機械判読可能なデータの表記方法の統一ルールの策定](https://www.soumu.go.jp/menu_news/s-news/01toukatsu01_02000186.html)というものを出している。
# これが最適な提案かはさておき、データの記述に法則性と機械判読性をもたせる意識を全員が持つことが重要なように思う。
# 
# お掃除ロボットが床を綺麗にするためには、まずお掃除ロボットが走れるよう掃除する(床に物が散乱していない)という条件が求められる、という話だ(そうなの？)。

# ## パスの指定
# 
# ファイルがコンピュータ上でどこにあるかを指し示す文字列はファイルパス(path)と呼ばれる。  
# ```"/content/drive/My Drive/XXX.png"```もファイルパスの一例になっている。
# 
# たとえば...  
# >[Sota]というユーザの[ドキュメント] (あるいは[書類])フォルダに  
# [csv_file]というフォルダがあり[test.csv]というcsvファイルが入っている
# 
# とするとそのファイルを指し示すパスは  
# Windowsの場合→ ```C:\Users\Sota\Douments\csv_file\test.csv```  
# macOSの場合→ ```/Users/Sota/Documents/csv_file/test.csv```
# となる。
# 
# 注:  
# * Windowsの場合→"C"の部分は皆さんのディスク環境に依存
# * Google Colab.環境では、Unix(Mac)やLinuxと同様の方式(スラッシュを使用)  
# * バックスラッシュ\はWindowsの日本語環境では¥円記号で表示される  
# (プログラムなどを書く人には厄介な仕様だったりする)  
# 
# コンピュータには、ホームディレクトリというものが指定されておりWindowsなら ```C:\Users\ユーザー名```,Macなら ```/Users/ユーザー名```に通常設定されていて、ユーザーがよく使うデスクトップや写真・ドキュメントなどのフォルダはホームディレクトリ直下に配置されている。また、ホームディレクトリは```~/```で簡略化して指定することもできる。
# OSにもよるが...ライトユーザーはホームディレクトリより上の階層をあまり触らないことが推奨されている(と思う)。理由は、システムファイルが入っていることが多いため。
# 
# パスの指定の仕方にはその他にも方法があり、ピリオドやスラッシュを駆使して現在のディレクトリからの[相対パス]で指定する事もできる。たとえば...
# 
# Home  
# ├ Documents  
# │└─ AdDS2021  
# ││   　└─ Report1  
# │└─ AdDS2020  
# ││   　└─ Report1  
# ││   　│　 └─ StudentA  
# ││   　│　 └─ StudentB  
# ││   　└─ Report2  
# │└─ AdDS2019  
# ├ Picures  
# ︙
# 
# こういう階層構造になっていて、現在```Home/Documents/AdDS2020/Report1```という
# ディレクトリにいるとすると、そこから
# * StudentAへの相対パスは ```./StudentA```
# * Report2への相対パスは ```../Report2```
# * AdDS2019への相対パスは ```../../AdDS2019```
# * Pictureへの相対パスは```../../../Pictures```
# 
# といった感じ。前述のように愚直にReport1フォルダを指定するときは```/Users/Sota/Documents/AdDS2020/Report1```といった感じで、これを相対パスと対比させて絶対パスと呼んだりする。

# ### 余談: ファイル名に使用すべきでない文字
# 
# 授業で公開しているノートブックの名前は基本的に半角英数字とアンダースコアのみで構成されている。これは別に作成者(吉田)がイキってる訳ではない。
# 
# * 半角スペース(以下␣という記号で表現する)
# * 各種括弧 (),{},[]
# * カンマ ,
# * ピリオド .
# * ハイフン -
# * スラッシュ /
# * エクスクラメーションマーク !
# * 円記号(バックスラッシュ) ¥
# * その他、機種依存文字などはもちろん、全角記号等
# 
# などは、(プログラムで扱う予定がある)ファイルの名前には使用しないことが推奨される。その理由は色々あるが
# 
# 1. 機械の解釈にambiguity(あいまいさ)が生じる
# 2. (1.により人間側の操作が増えて)面倒
# 
# というところに尽きる。例を示そう。  
# Google Colab.上では冒頭に!を付けることで、以下に例を示すようなLinuxコマンドを実行できる。
# 
# ```!ls hogehoge.pdf``` #← lsコマンド　リスト(該当ファイル等)を表示  
# ```!mkdir hogehoge``` #← make directoryコマンド  
# ```!rm hogehoge``` #←remove(削除)コマンド  

# たとえば半角スペースが入った```test␣.pdf```というファイルがあったとする。  
# これをlsコマンドで表示させようとして
# 
# ```
# !ls test .pdf
# ```
# 
# という命令を行うと、```test␣.pdf```という指定ではなく```test```と```.pdf```という2つのファイルが存在するかどうかを尋ねる命令になってしまう。  
# この場合、```test␣.pdf```の有無を調べたければ、別途バックスラッシュを入れて「記号としての空白です」と機械に教えなくてはならない。  
# ```
# !ls test\ .pdf
# ```
# 
# といった具合に、人間側の手間が必要になってしまう。  
# 人間が目で見るフォルダ名と機械に与えるべきパスが異なるというのは...やっぱり色んな場面で不便が生じる。
# 上記の記号や２バイト以上の文字はファイル(フォルダ)名に使わないのがコンピューターにとっては無難なのだ。
# こういうことは小中高や大学でも理由付きで教えてくれなかったりするので、プログラミングをやって初めて気がつく(気にするようになった)という人も多いかもしれない。
# 
