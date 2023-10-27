p1 = [2.0, 4.0, -5.0]
p2 = [1.0, 3.0, -4.0]
#...中略
p100 = [5.5,-2.0, 3.0]

d_1_2 = ( (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2 ) ** 0.5
d_1_100 = ( (p1[0] - p100[0])**2 + (p1[1] - p100[1])**2 + (p1[2] - p100[2])**2 ) ** 0.5

def calc_d(l1,l2):           
    return ( (l1[0] - l2[0])**2 + (l1[1] - l2[1])**2 + (l1[2] - l2[2])**2 ) ** 0.5

t = calc_d(p1,p2) 
print("点1",p1, "と点2", p2, "の距離は", t, "です")

print(calc_d(p1,p100)) #←これでも使えるし
print(calc_d([20.0, 1.0,-5.0], [-2.0, 3.0,5.5])) #←などとして名前をつけなくても使える

import random 
lists = [ [ random.gauss(0,1) for n in range(3)] for i in range(100)] #3次元の座標点をランダムに100個作っている n,iはダミー変数(特に使ってない)
hit = 0
for j in range(100):
    for i in range(j+1,100): # i>j
        distance = calc_d( lists[j], lists[i])
        #print(j,i, distance) # 4950回文の計算結果をprintすると邪魔なのでコメントアウトした
        hit += 1 
print(hit) #回数だけ表示しよう
#上のjのループ内で、iはj+1から99までを回る。 j+1= 100つまり j=99のとき range(j+1,100)はちゃんと空になる
#つまり、長さ100のリストにindex=100でアクセス(範囲外参照)したりすることはない。

def name(): #引数なしで、ただ以下の文字列を表示する関数
    print("私は田中です")

def myname(namae): #引数namaeを使って、以下の文字列を表示する関数
    print("私は"+str(namae)+"です")

def myname_return(namae): #　myname()で表示させた文字列自体を返す関数
    return "私は"+str(namae)+"です"

print("name()の実行→", name()) ## name()が実行されたあとにココのprint文が実行される。

print("myname()の返り値→", myname("吉田"))

print("myname_return()の返り値→", myname_return("吉田"))

def calc_d_print(l1,l2):
    return "距離は"+str( ( (l1[0] - l2[0])**2 + (l1[1] - l2[1])**2 + (l1[2] - l2[2])**2 ) ** 0.5  )+"です"

def zahyo_and_d(l1,l2):
    d = calc_d(l1,l2) #関数の中で、先程の自作関数を呼んでいる
    return [l1,l2],d  #座標を結合したリストと距離を返す

ret = calc_d_print(p1,p2)
print("関数calc_d_print→", ret,type(ret))


ret = zahyo_and_d(p1,p2)
print("関数zahyo→ ", ret,type(ret))
print("座標の結合リスト",ret[0],"距離",ret[1])

#数値リストの要素のp乗和を計算する関数
def sump(tmp,p=2): 
    return sum([tmp[i]**p for i in range(len(tmp))])

list1 = [10.0,20.0,30.0,40.0]
print("default", sump(list1)) #pを指定しなければp=2が選ばれる
print("p=1", sump(list1,p=1))
print("p=2", sump(list1,2))
print("p=3", sump(list1,3))

a = 2
list1 = [10.0,20.0,30.0,40.0]

def testfunc():
    print(a)

a = 2
testfunc()

def testfunc():
    abcd = 1.000
testfunc()
print(abcd)

def testfunc():
    a = 5
a= 2
testfunc()
print(a)

def testfunc():
    a = 5
    print("関数の内部", a, id(a))
    
a= 2 
print("関数の実行前", a, id(a))
testfunc()
print("関数の実行後", a, id(a)) 

def testfunc():
    global abc, a #global変数の宣言
    abc = 5
    a += 2

a=2
print("実行前")
print("a",a , id(a))
testfunc()
print("実行後")
print("a", a, id(a))  #別の変数として再定義されていることが分かる
print("abc", abc)

def func_join(listA,listB): #特殊なリストの結合をして返す関数
    return listA + 2 * listB 

list1 = [ 2.0,30.0,18.0]
list2 = [ 9.0,4.0,8.0]
nlist = func_join(list1,list2)

def func_update_list(in_list):
    in_list[0] = "AAA"

tmp = [ "SS", 1,2,3]
print("実行前", tmp, id(tmp), id(tmp[0]))
func_update_list(tmp)
print("実行後", tmp,id(tmp),id(tmp[0])) 

def func_update_list(in_list):
    in_list[0] = "AAA" 
    in_list = ["BBB", 0,1,2]  ##ココはローカル変数扱い
    return in_list

tmp = [ "SS", 1,2,3]
print("実行前", tmp, id(tmp), id(tmp[0]))
ret = func_update_list(tmp)
print("実行後", tmp,id(tmp),id(tmp[0])) 
print("ret", ret,id(ret),id(ret[0])) 

a = [1,2,3]
a.append(4)
print(a)
